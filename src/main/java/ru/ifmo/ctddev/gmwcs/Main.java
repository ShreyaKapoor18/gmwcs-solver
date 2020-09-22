package ru.ifmo.ctddev.gmwcs;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import ru.ifmo.ctddev.gmwcs.graph.*;
import ru.ifmo.ctddev.gmwcs.solver.BicomponentSolver;
import ru.ifmo.ctddev.gmwcs.solver.RLTSolver;
import ru.ifmo.ctddev.gmwcs.solver.Solver;
import ru.ifmo.ctddev.gmwcs.solver.SolverException;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.List;

import static java.util.Arrays.asList;

public class Main {
    public static OptionSet parseArgs(String args[]) throws IOException {
        OptionParser optionParser = new OptionParser();
        optionParser.allowsUnrecognizedOptions();
        optionParser.acceptsAll(asList("h", "help"), "Print a short help message");
        optionParser.accepts("version");
        OptionSet optionSet = optionParser.parse(args);
        optionParser.acceptsAll(asList("n", "nodes"), "Node list file").withRequiredArg().required();
        optionParser.acceptsAll(asList("e", "edges"), "Edge list file").withRequiredArg().required();
        optionParser.acceptsAll(asList("nm", "number"), "Preserved Nodes").withRequiredArg()
                .ofType(Double.class).defaultsTo(10.0);
        optionParser.accepts("root", "Root node").withRequiredArg();
        optionParser.acceptsAll(asList("m", "threads"), "Number of threads").withRequiredArg()
                .ofType(Integer.class).defaultsTo(1);
        optionParser.acceptsAll(asList("t", "timelimit"), "Timelimit in seconds (<= 0 - unlimited)")
                .withRequiredArg().ofType(Long.class).defaultsTo(0L);
        optionParser.acceptsAll(asList("u", "unrooted"), "Maximum share of time allocated for solving unrooted parts")
                .withRequiredArg().ofType(Double.class).defaultsTo(0.3);
        optionParser.acceptsAll(asList("r", "rooted"), "Maximum share of time allocated for solving rooted parts")
                .withRequiredArg().ofType(Double.class).defaultsTo(0.3);
        // add option for the max number of nodes here
        if (optionSet.has("h")) {
            optionParser.printHelpOn(System.out);
            return null;
        }
        if (optionSet.has("version")) {
            System.out.println("gmwcs-solver version " + Main.class.getPackage().getImplementationVersion());
            return null;
        }
        try {
            optionSet = optionParser.parse(args);
            double ush = (Double) optionSet.valueOf("u");
            double rsh = (Double) optionSet.valueOf("r");
            if (ush < 0.0 || ush > 1.0 || rsh < 0.0 || rsh > 1.0) {
                System.err.println("Share must b in range [0,1]");
                return null;
            }
            if (rsh + ush > 1.0) {
                System.err.println("Sum of shares of rooted and unrooted parts must be <= 1.0");
                return null;
            }
        } catch (Exception e) {
            System.err.println(e.getMessage());
            System.err.println();
            optionParser.printHelpOn(System.err);
            return null;
        }
        return optionSet;
    }

    public static void main(String[] args) {
        OptionSet optionSet = null;
        try {
            optionSet = parseArgs(args);
            if (optionSet == null) {
                return;
            }
        } catch (IOException e) {
            // We can't say anything. Error occurred while printing to stderr.
            return;
        }
        long timelimit = (Long) optionSet.valueOf("timelimit");
        TimeLimit tl = new TimeLimit(timelimit <= 0 ? Double.POSITIVE_INFINITY : timelimit);
        double rsh = (Double) optionSet.valueOf("r");
        double ush = (Double) optionSet.valueOf("u");
        double max_num_nodes = (Double) optionSet.valueOf("nm");
        TimeLimit biggestTL = tl.subLimit(1.0 - ush);
        int threadsNum = (Integer) optionSet.valueOf("threads");
        File nodeFile = new File((String) optionSet.valueOf("nodes"));
        File edgeFile = new File((String) optionSet.valueOf("edges"));
        RLTSolver rltSolver = new RLTSolver(max_num_nodes); // first entry point for the solver
        rltSolver.setThreadsNum(threadsNum);
        Solver solver = null;
        if (!optionSet.has("root")) {
            BicomponentSolver comp_solver = new BicomponentSolver(rltSolver);
            comp_solver.setUnrootedTL(tl);
            comp_solver.setRootedTL(biggestTL.subLimit(ush == 1.0 ? 0 : rsh / (1.0 - ush)));
            comp_solver.setTLForBiggest(biggestTL);
            solver = comp_solver;
        } else {
            solver = rltSolver;
            solver.setTimeLimit(tl);
        }
        GraphIO graphIO = new SimpleIO(nodeFile, new File(nodeFile.toString() + ".out"),
                edgeFile, new File(edgeFile.toString() + ".out"));
        try {
            Graph graph = graphIO.read();
            if (optionSet.has("root")) {
                Node root = graphIO.nodeByName((String) optionSet.valueOf("root"));
                if (root == null) {
                    System.err.println("Chosen root node is not presented in the graph");
                    return;
                }
                rltSolver.setRoot(root);
            }
            List<Unit> units = solver.solve(graph); // if no root options then BiComponent, if rooted then rltsolver
            graphIO.write(units);
        } catch (ParseException e) {
            System.err.println("Couldn't parse input files: " + e.getMessage() + " " + e.getErrorOffset());
        } catch (SolverException e) {
            System.err.println("Error occur while solving:" + e.getMessage());
        } catch (IOException e) {
            System.err.println("Error occurred while reading/writing input/output files");
        }
    }
}
