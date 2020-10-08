package ru.ifmo.ctddev.gmwcs.solver;

import ru.ifmo.ctddev.gmwcs.Pair;
import ru.ifmo.ctddev.gmwcs.TimeLimit;
import ru.ifmo.ctddev.gmwcs.graph.*;

import java.util.*;

public class BicomponentSolver implements Solver { //bicomponent solver is the child of Solver
    private TimeLimit rooted;
    private TimeLimit biggest;
    private TimeLimit unrooted;
    private RLTSolver solver;
    private boolean isSolvedToOptimality;
    private double lb;
    private boolean silence;

    public BicomponentSolver(RLTSolver solver) {
        //rooted = new TimeLimit(Double.POSITIVE_INFINITY);
        rooted = new TimeLimit(10000);
        unrooted = biggest = rooted;
        this.solver = solver;
        lb = 0;
    }

    public void setRootedTL(TimeLimit tl) {
        this.rooted = tl;
    }

    public void setUnrootedTL(TimeLimit tl) {
        this.unrooted = tl;
    }

    public void setTLForBiggest(TimeLimit tl) {
        this.biggest = tl;
    }

    public List<Unit> solve(Graph graph) throws SolverException {
        // in the unrooted case this will be executed
        Graph g = graph;
        graph = graph.subgraph(graph.vertexSet()); //how is it extracting the subpgraph
        //Preprocessor.preprocess(graph); //silencing the preprocessor
        //if (!silence) {
        //    System.out.print("Preprocessing deleted " + (g.vertexSet().size() - graph.vertexSet().size()) + " nodes ");
         //   System.out.println("and " + (g.edgeSet().size() - graph.edgeSet().size()) + " edges.");
        //}
        System.out.println("solve method of bicomponent solver");
        isSolvedToOptimality = true;
        solver.setLB(-Double.MAX_VALUE);
        if (graph.vertexSet().size() == 0) {
            return null;
        }
        long timeBefore = System.currentTimeMillis();
        Decomposition decomposition = new Decomposition(graph);
        double duration = (System.currentTimeMillis() - timeBefore) / 1000.0;
        if (!silence) {
            System.out.println("Graph decomposing takes " + duration + " seconds.");
        }
        System.out.println("decomposed part of the graph" + decomposition);
        List<Unit> bestBiggest = solveBiggest(graph, decomposition);
        List<Unit> bestUnrooted = extract(solveUnrooted(graph, decomposition));
        graph.vertexSet().forEach(Node::clear);
        graph.edgeSet().forEach(Edge::clear);
        List<Unit> best = Utils.sum(bestBiggest) > Utils.sum(bestUnrooted) ? bestBiggest : bestUnrooted;
        solver.setLB(-Double.MAX_VALUE);
        if (Utils.sum(best) < 0) {
            return null;
        }
        return best;
    }

    private List<Unit> extract(List<Unit> sol) {
        System.out.println("extract method of bicomponent solver");
        List<Unit> res = new ArrayList<>();
        for (Unit u : sol) {
            for (Unit a : u.getAbsorbed()) {
                res.add(a);
            }
            res.add(u);
        }
        return res;
    }

    @Override
    public void setTimeLimit(TimeLimit tl) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isSolvedToOptimality() {
        return isSolvedToOptimality;
    }

    @Override
    public void suppressOutput() {
        solver.suppressOutput();
        silence = true;
    }

    @Override
    public void setLB(double lb) {
        this.lb = lb;
    }

    private Node getRoot(Graph graph) {
        Set<Node> rootCandidates = new LinkedHashSet<>();
        for (int i = -1; i < graph.vertexSet().size(); i++) {
            rootCandidates.add(new Node(i, 0.0)); // here it is looking for the root vertex
        }
        graph.vertexSet().stream().forEach(v -> rootCandidates.removeAll(v.getAbsorbed())); // but we are not getting any absorbed vertices
        rootCandidates.removeAll(graph.vertexSet());
        return rootCandidates.iterator().next();
    }

    private List<Unit> solveBiggest(Graph graph, Decomposition decomposition) throws SolverException {
        System.out.println("solve biggest");
        Graph tree = new Graph();
        Map<Unit, List<Unit>> history = new HashMap<>();
        graph.vertexSet().forEach(v -> history.put(v, v.getAbsorbed())); // but the preprocessing is disabled
        graph.edgeSet().forEach(e -> history.put(e, e.getAbsorbed())); // need to change this since the processing
        Node root = getRoot(graph);
        tree.addVertex(root); // first add only the root vertex for the -1
        System.out.println("root vertex added"+ root);
        Map<Unit, Node> itsCutpoints = new LinkedHashMap<>();
        for (Pair<Set<Node>, Node> p : decomposition.getRootedComponents()) {
            for (Node node : p.first) {
                for (Edge edge : graph.edgesOf(node)) {
                    itsCutpoints.put(edge, p.second);
                }
                itsCutpoints.put(node, p.second);
            }
            tree.addGraph(graph.subgraph(p.first));
            addAsChild(tree, p.first, p.second, root);
        }
        solver.setRoot(root);
        System.out.println("trying to solve tree with the rooted part");
        List<Unit> rootedRes = solve(tree, rooted);
        solver.setRoot(null);
        Graph main = graph.subgraph(decomposition.getBiggestComponent());
        if (rootedRes != null) {
            rootedRes.stream().filter(unit -> unit != root).forEach(unit -> {
                Node cutpoint = itsCutpoints.get(unit);
                cutpoint.absorb(unit);
            });
        }
        solver.setLB(lb);
        System.out.println("solve for the main vs biggest");
        System.out.println("main:" + main);
        System.out.println("biggest"+ biggest);
        List<Unit> solution = solve(main, biggest);
        List<Unit> result = new ArrayList<>();
        result.addAll(solution);
        solver.setLB(Utils.sum(result));
        solution.stream().forEach(u -> result.addAll(u.getAbsorbed()));
        repairCutpoints(history);
        return result;
    }

    private void repairCutpoints(Map<Unit, List<Unit>> history) {
        history.keySet().forEach(Unit::clear);
        for (Unit u : history.keySet()) {
            for (Unit a : history.get(u)) {
                u.absorb(a);
            }
        }
    }

    private void addAsChild(Graph tree, Set<Node> component, Node cp, Node root) {
        for (Node neighbour : tree.neighborListOf(cp)) {
            if (!component.contains(neighbour)) {
                continue;
            }
            Edge edge = tree.getEdge(cp, neighbour);
            tree.removeEdge(edge);
            tree.addEdge(root, neighbour, edge);
        }
        tree.removeVertex(cp);
    }

    private List<Unit> solveUnrooted(Graph graph, Decomposition decomposition) throws SolverException {
        Set<Node> union = new LinkedHashSet<>();
        decomposition.getUnrootedComponents().forEach(union::addAll);
        return solve(graph.subgraph(union), unrooted);
    }

    private List<Unit> solve(Graph graph, TimeLimit tl) throws SolverException {
        System.out.println("solving with the time limit on the graph");
        solver.setTimeLimit(tl);
        System.out.println("Trying to use the solve method");
        List<Unit> result = solver.solve(graph);
        if (!solver.isSolvedToOptimality()) {
            isSolvedToOptimality = false;
        }
        return result;
    }
}
