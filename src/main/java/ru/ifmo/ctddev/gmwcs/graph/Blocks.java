package ru.ifmo.ctddev.gmwcs.graph;

import org.jgrapht.Graphs;
import org.jgrapht.UndirectedGraph;

import java.util.*;

public class Blocks {
    private UndirectedGraph<Node, Edge> graph;
    private Map<Node, Integer> enter;
    private Map<Node, Integer> up;
    private Stack<Edge> stack;
    private Set<Set<Node>> components;
    private Set<Node> cutpoints;
    private Node root;
    private int rootChildren;
    private int time;

    public Blocks(UndirectedGraph<Node, Edge> graph) {
        enter = new LinkedHashMap<>();
        up = new LinkedHashMap<>();
        stack = new Stack<>();
        root = graph.vertexSet().iterator().next();
        components = new LinkedHashSet<>();
        cutpoints = new LinkedHashSet<>();
        this.graph = graph;
        if (graph.vertexSet().size() > 1) {
            dfs(root, null);
        } else {
            Set<Node> component = new LinkedHashSet<>();
            component.add(root);
            components.add(component);
        }
    }

    public void dfs(Node v, Node parent) {
        time++;
        enter.put(v, time);
        up.put(v, time);
        for (Node u : Graphs.neighborListOf(graph, v)) {
            if (u == parent) {
                continue;
            }
            if (!enter.containsKey(u)) {
                stack.add(graph.getEdge(v, u));
                if (v == root) {
                    ++rootChildren;
                }
                dfs(u, v);
                if (up.get(u) >= enter.get(v)) {
                    Set<Node> component = new LinkedHashSet<>();
                    Edge expected = graph.getEdge(v, u);
                    while (true) {
                        Edge edge = stack.pop();
                        component.add(graph.getEdgeSource(edge));
                        component.add(graph.getEdgeTarget(edge));
                        if (edge == expected) {
                            break;
                        }
                    }
                    components.add(component);
                    cutpoints.add(v);
                }
                if (up.get(u) < up.get(v)) {
                    up.put(v, up.get(u));
                }
            } else {
                if (up.get(v) > enter.get(u)) {
                    up.put(v, enter.get(u));
                }
            }
        }
        if (rootChildren < 2) {
            cutpoints.remove(root);
        }
    }

    public Set<Set<Node>> components() {
        return components;
    }

    public Set<Node> cutpoints() {
        return cutpoints;
    }
}