import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.geom.Line2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;
import java.util.List;


public class Smart_Map_Navigation {
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> new SmartMapFrame(new InMemoryBackendService()).setVisible(true));
    }
}


interface BackendService {
    Graph loadGraph();
    void updateEdgeFactors(String from, String to, double traffic, double risk, boolean blocked);
}

class InMemoryBackendService implements BackendService {
    private final Graph graph;

    public InMemoryBackendService() {
        this.graph = DemoGraphs.cityDemo();
    }

    @Override
    public Graph loadGraph() {
        return graph;
    }

    @Override
    public void updateEdgeFactors(String from, String to, double traffic, double risk, boolean blocked) {
        Edge e = graph.getEdge(from, to);
        if (e != null) {
            e.traffic = traffic;
            e.risk = risk;
            e.blocked = blocked;
        }
    }
}


class Graph {
    final Map<String, Node> nodes = new LinkedHashMap<>();
    final Map<String, List<Edge>> adj = new LinkedHashMap<>();

    public void addNode(Node n) {
        nodes.put(n.id, n);
        adj.computeIfAbsent(n.id, k -> new ArrayList<>());
    }

    public void addEdge(String from, String to, double distanceKm, double travelTimeMin, double traffic, double risk, boolean blocked) {
        Node a = nodes.get(from); Node b = nodes.get(to);
        if (a == null || b == null) throw new IllegalArgumentException("Unknown node in edge");
        Edge e1 = new Edge(a, b, distanceKm, travelTimeMin, traffic, risk, blocked);
        Edge e2 = new Edge(b, a, distanceKm, travelTimeMin, traffic, risk, blocked); // undirected
        adj.get(from).add(e1);
        adj.get(to).add(e2);
    }

    public List<Edge> neighbors(String id) { return adj.getOrDefault(id, Collections.emptyList()); }

    public Edge getEdge(String from, String to) {
        List<Edge> list = adj.get(from);
        if (list == null) return null;
        for (Edge e : list) if (e.to.id.equals(to)) return e;
        return null;
    }
}

class Node {
    final String id;
    final String label;
    final int x, y; 

    public Node(String id, String label, int x, int y) {
        this.id = id; this.label = label; this.x = x; this.y = y;
    }

    @Override public String toString() { return label; }
}

class Edge {
    final Node from, to;
    final double distanceKm;         
    final double travelTimeMin;      
    double traffic;                  
    double risk;                     
    boolean blocked;                 

    public Edge(Node from, Node to, double distanceKm, double travelTimeMin, double traffic, double risk, boolean blocked) {
        this.from = from; this.to = to;
        this.distanceKm = distanceKm;
        this.travelTimeMin = travelTimeMin;
        this.traffic = clamp01(traffic);
        this.risk = clamp01(risk);
        this.blocked = blocked;
    }

    static double clamp01(double v) { return Math.max(0, Math.min(1, v)); }
}


class Dijkstra {
    static class Result {
        final Map<String, Double> dist;
        final Map<String, String> prev;
        Result(Map<String, Double> d, Map<String, String> p) { this.dist = d; this.prev = p; }
        List<String> path(String target) {
            LinkedList<String> path = new LinkedList<>();
            String cur = target; if (!prev.containsKey(cur) && !dist.containsKey(cur)) return path;
            while (cur != null) { path.addFirst(cur); cur = prev.get(cur); }
            return path;
        }
    }

    public static Result shortestPath(Graph g, String source, WeightModel wm) {
        Map<String, Double> dist = new HashMap<>();
        Map<String, String> prev = new HashMap<>();
        Set<String> visited = new HashSet<>();
        PriorityQueue<String> pq = new PriorityQueue<>(Comparator.comparingDouble(dist::get));

        for (String id : g.nodes.keySet()) dist.put(id, Double.POSITIVE_INFINITY);
        dist.put(source, 0.0);
        pq.add(source);

        while (!pq.isEmpty()) {
            String u = pq.poll();
            if (visited.contains(u)) continue;
            visited.add(u);

            for (Edge e : g.neighbors(u)) {
                if (e.blocked) continue; // skip closed edges
                double w = wm.weight(e);
                if (w == Double.POSITIVE_INFINITY) continue;
                double alt = dist.get(u) + w;
                if (alt < dist.get(e.to.id)) {
                    dist.put(e.to.id, alt);
                    prev.put(e.to.id, u);
                    pq.add(e.to.id);
                }
            }
        }
        return new Result(dist, prev);
    }
}


interface WeightModel {
    double weight(Edge e);
}

class CompositeWeight implements WeightModel {
    private final double alpha; 
    private final double beta; 
    private final double gamma; 
    private final double delta; 
    private final double blockPenalty;

    public CompositeWeight(double alpha, double beta, double gamma, double delta, double blockPenalty) {
        this.alpha = alpha; this.beta = beta; this.gamma = gamma; this.delta = delta; this.blockPenalty = blockPenalty;
    }

    @Override
    public double weight(Edge e) {
        if (e.blocked) return Double.POSITIVE_INFINITY;
        double trafficFactor = 1.0 + 2.0 * e.traffic; 
        double timeAdj = e.travelTimeMin * trafficFactor;
        double nearBlock = (e.risk > 0.85 || e.traffic > 0.9) ? blockPenalty : 0.0;
        return alpha * e.distanceKm + beta * timeAdj + gamma * (e.risk * 100.0) + delta * (e.traffic * 100.0) + nearBlock;
    }
}


class DemoGraphs {
    static Graph cityDemo() {
        Graph g = new Graph();
        // Nodes laid out like a small city map (800x520 canvas)
        g.addNode(new Node("A", "A — Mumbai", 90, 420));
        g.addNode(new Node("B", "B — Delhi", 220, 440));
        g.addNode(new Node("C", "C — Bangalore", 360, 460));
        g.addNode(new Node("D", "D —Kolkata", 520, 440));
        g.addNode(new Node("E", "E —  Chennai", 660, 420));
        g.addNode(new Node("F", "F —  Hyderabad", 520, 300));
        g.addNode(new Node("G", "G — Ahmedabad", 360, 240));
        g.addNode(new Node("H", "H — Pune", 220, 280));
        g.addNode(new Node("S", "S — Surat", 660, 240));

       
        g.addEdge("A","B", 2.0, 6, 0.2, 0.1, false);
        g.addEdge("B","C", 1.6, 5, 0.3, 0.2, false);
        g.addEdge("C","D", 1.8, 5, 0.5, 0.3, false);
        g.addEdge("D","E", 2.0, 6, 0.2, 0.1, false);
        g.addEdge("B","H", 2.2, 7, 0.4, 0.1, false);
        g.addEdge("H","G", 2.0, 6, 0.5, 0.2, false);
        g.addEdge("G","F", 2.0, 6, 0.3, 0.2, false);
        g.addEdge("F","D", 2.0, 6, 0.2, 0.2, false);
        g.addEdge("C","G", 2.2, 8, 0.6, 0.4, false);
        g.addEdge("F","S", 2.0, 6, 0.2, 0.1, false);
        g.addEdge("D","S", 2.4, 7, 0.3, 0.2, false);
        g.addEdge("E","S", 2.2, 6, 0.2, 0.1, false);
        g.addEdge("A","H", 2.6, 9, 0.3, 0.2, false);
        g.addEdge("B","G", 2.8, 10, 0.6, 0.3, false);
        // Potentially risky riverside segment
        g.addEdge("C","E", 3.4, 12, 0.4, 0.7, false);
        return g;
    }
}

/* ********************************* UI ********************************/
class SmartMapFrame extends JFrame {
    private final BackendService backend;
    private final Graph graph;

    private final JComboBox<Node> startBox;
    private final JComboBox<Node> endBox;
    private final JButton routeBtn;
    private final JButton randomizeBtn;
    private final JCheckBox avoidClosedBox;
    private final JSlider alphaSlider, betaSlider, gammaSlider, deltaSlider;
    private final JLabel statusLabel;
    private final JTextArea infoArea;
    private final MapPanel mapPanel;

    // For drawing selected path
    private List<String> currentPath = new ArrayList<>();

    public SmartMapFrame(BackendService backend) {
        super("Smart Map Navigation — Dijkstra");
        this.backend = backend;
        this.graph = backend.loadGraph();

        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setSize(1024, 720);
        setLocationRelativeTo(null);
        setLayout(new BorderLayout());

        // Left control panel
        JPanel left = new JPanel(new BorderLayout());
        left.setPreferredSize(new Dimension(320, 720));
        left.setBorder(new EmptyBorder(12, 12, 12, 12));
        add(left, BorderLayout.WEST);

        // Top: selectors
        JPanel selectors = new JPanel();
        selectors.setLayout(new GridLayout(0,1,8,8));
        startBox = new JComboBox<>(graph.nodes.values().toArray(new Node[0]));
        endBox = new JComboBox<>(graph.nodes.values().toArray(new Node[0]));
        routeBtn = new JButton("Find Route");
        randomizeBtn = new JButton("Simulate Disaster/Traffic");
        avoidClosedBox = new JCheckBox("Auto-close very risky roads (risk>0.85 or traffic>0.9)", true);

        selectors.add(new JLabel("Start:")); selectors.add(startBox);
        selectors.add(new JLabel("Destination:")); selectors.add(endBox);
        selectors.add(routeBtn); selectors.add(randomizeBtn); selectors.add(avoidClosedBox);

        left.add(selectors, BorderLayout.NORTH);

        // Middle: weight sliders
        JPanel sliders = new JPanel();
        sliders.setLayout(new GridLayout(0,1,6,6));
        alphaSlider = labeledSlider(sliders, "α Distance weight", 0, 100, 20);
        betaSlider  = labeledSlider(sliders, "β Time weight",     0, 100, 40);
        gammaSlider = labeledSlider(sliders, "γ Risk weight",     0, 100, 50);
        deltaSlider = labeledSlider(sliders, "δ Traffic weight",  0, 100, 30);
        left.add(sliders, BorderLayout.CENTER);

        // Bottom
        JPanel bottom = new JPanel(new BorderLayout());
        statusLabel = new JLabel("Ready.");
        infoArea = new JTextArea(8, 30);
        infoArea.setEditable(false);
        infoArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        JScrollPane sp = new JScrollPane(infoArea);
        bottom.add(statusLabel, BorderLayout.NORTH);
        bottom.add(sp, BorderLayout.CENTER);
        left.add(bottom, BorderLayout.SOUTH);

        // Map panel
        mapPanel = new MapPanel(graph, () -> currentPath);
        add(mapPanel, BorderLayout.CENTER);

        // Actions
        routeBtn.addActionListener(new RouteAction());
        randomizeBtn.addActionListener(e -> simulateDisaster());

        ChangeListener live = new ChangeListener() {
            @Override public void stateChanged(ChangeEvent e) { findRoute(false); }
        };
        alphaSlider.addChangeListener(live);
        betaSlider.addChangeListener(live);
        gammaSlider.addChangeListener(live);
        deltaSlider.addChangeListener(live);

        // Initial path
        startBox.setSelectedItem(graph.nodes.get("A"));
        endBox.setSelectedItem(graph.nodes.get("S"));
        findRoute(true);
    }

    private JSlider labeledSlider(JPanel host, String label, int min, int max, int val) {
        JPanel p = new JPanel(new BorderLayout());
        JLabel l = new JLabel(label + " (" + val + ")");
        JSlider s = new JSlider(min, max, val);
        s.addChangeListener(e -> l.setText(label + " (" + s.getValue() + ")"));
        p.add(l, BorderLayout.NORTH);
        p.add(s, BorderLayout.CENTER);
        host.add(p);
        return s;
    }

    private void simulateDisaster() {
        Random r = new Random();
        for (String from : graph.adj.keySet()) {
            for (Edge e : graph.neighbors(from)) {
                // Only update forward once per pair to avoid double-random in undirected edges
                if (e.from.id.compareTo(e.to.id) < 0) {
                    double newTraffic = Math.min(1.0, Math.max(0.0, e.traffic + (r.nextGaussian()*0.2)));
                    double newRisk = Math.min(1.0, Math.max(0.0, e.risk + (r.nextGaussian()*0.2)));
                    boolean blocked = avoidClosedBox.isSelected() && (newRisk > 0.88 || newTraffic > 0.92);
                    backend.updateEdgeFactors(e.from.id, e.to.id, newTraffic, newRisk, blocked);
                }
            }
        }
        infoArea.append("\n[SIM] Updated edges with random traffic/risk; near-critical segments may be closed.\n");
        findRoute(true);
    }

    private void findRoute(boolean log) {
        Node start = (Node) startBox.getSelectedItem();
        Node end = (Node) endBox.getSelectedItem();
        if (start == null || end == null) return;

        double alpha = alphaSlider.getValue() / 10.0;
        double beta  = betaSlider.getValue()  / 10.0;
        double gamma = gammaSlider.getValue() / 10.0;
        double delta = deltaSlider.getValue() / 10.0;

        CompositeWeight wm = new CompositeWeight(alpha, beta, gamma, delta, 500.0);
        Dijkstra.Result res = Dijkstra.shortestPath(graph, start.id, wm);
        List<String> path = res.path(end.id);
        currentPath = path;

        if (path.isEmpty()) {
            statusLabel.setText("No available path (some roads may be blocked).");
            infoArea.append("No path from " + start.label + " to " + end.label + " with current conditions.\n");
        } else {
            double cost = res.dist.get(end.id);
            statusLabel.setText("Path found: cost=" + String.format(Locale.US, "%.2f", cost) + ", hops=" + (path.size()-1));
            if (log) {
                infoArea.append("\nRoute " + start.id + " → " + end.id + ": " + path + "\n");
                infoArea.append(describePath(path));
            }
        }
        mapPanel.repaint();
    }

    private String describePath(List<String> path) {
        StringBuilder sb = new StringBuilder();
        double totalDist = 0, totalTime = 0; double maxRisk = 0, maxTraffic = 0; boolean anyBlocked = false;
        for (int i=0;i<path.size()-1;i++) {
            Edge e = graph.getEdge(path.get(i), path.get(i+1));
            if (e == null) continue;
            totalDist += e.distanceKm;
            double t = e.travelTimeMin * (1 + 2*e.traffic);
            totalTime += t;
            maxRisk = Math.max(maxRisk, e.risk);
            maxTraffic = Math.max(maxTraffic, e.traffic);
            anyBlocked |= e.blocked;
            sb.append(String.format(Locale.US,
                    "  %s → %s  dist=%.1f km  time≈%.1f min  traffic=%.2f  risk=%.2f%s\n",
                    e.from.id, e.to.id, e.distanceKm, t, e.traffic, e.risk, e.blocked?" [BLOCKED]":""));
        }
        sb.append(String.format(Locale.US,
                "TOTAL: distance=%.1f km  ETA≈%.1f min  peakTraffic=%.2f  peakRisk=%.2f%s\n",
                totalDist, totalTime, maxTraffic, maxRisk, anyBlocked?"  (contains blocked!)":""));
        return sb.toString();
    }

    private class RouteAction implements ActionListener {
        @Override public void actionPerformed(ActionEvent e) { findRoute(true); }
    }
}

class MapPanel extends JPanel {
    private final Graph graph;
    private final PathSupplier supplier;
    private final Stroke normalStroke = new BasicStroke(2f);
    private final Stroke pathStroke = new BasicStroke(5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

    public MapPanel(Graph graph, PathSupplier supplier) {
        this.graph = graph;
        this.supplier = supplier;
        setBackground(new Color(245, 247, 250));
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g.create();
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        // draw edges
        Set<String> pathPairs = new HashSet<>();
        List<String> path = supplier.currentPath();
        for (int i=0;i<path.size()-1;i++) pathPairs.add(pair(path.get(i), path.get(i+1)));

        for (String from : graph.adj.keySet()) {
            for (Edge e : graph.neighbors(from)) {
                if (!e.from.id.equals(from)) continue; // draw each undirected edge once
                if (e.from.id.compareTo(e.to.id) > 0) continue;

                Point p1 = new Point(e.from.x, e.from.y);
                Point p2 = new Point(e.to.x, e.to.y);

                // Edge color by status
                Color col = e.blocked ? new Color(180, 60, 60) : new Color(150, 160, 180);
                if (pathPairs.contains(pair(e.from.id, e.to.id))) col = new Color(60,120,200);

                g2.setStroke(normalStroke);
                g2.setColor(col);
                g2.draw(new Line2D.Double(p1, p2));

                // small label with traffic/risk
                String lbl = String.format(Locale.US, "t=%.2f r=%.2f", e.traffic, e.risk);
                int mx = (p1.x + p2.x)/2, my = (p1.y + p2.y)/2;
                g2.setFont(getFont().deriveFont(Font.PLAIN, 11f));
                g2.setColor(new Color(20,20,20,160));
                g2.drawString(lbl, mx+6, my-6);
            }
        }

        // highlight path on top
        g2.setStroke(pathStroke);
        g2.setColor(new Color(60,120,200));
        for (int i=0;i<path.size()-1;i++) {
            Node a = graph.nodes.get(path.get(i));
            Node b = graph.nodes.get(path.get(i+1));
            g2.draw(new Line2D.Double(a.x, a.y, b.x, b.y));
        }

        // draw nodes
        for (Node n : graph.nodes.values()) {
            int r = 18;
            int x = n.x - r/2, y = n.y - r/2;
            g2.setColor(new Color(255,255,255));
            g2.fillOval(x, y, r, r);
            g2.setColor(new Color(100, 110, 130));
            g2.setStroke(new BasicStroke(2f));
            g2.drawOval(x, y, r, r);
            g2.setColor(new Color(30,30,30));
            g2.setFont(getFont().deriveFont(Font.BOLD, 12f));
            g2.drawString(n.id, n.x - 4, n.y - 22);
            g2.setFont(getFont().deriveFont(Font.PLAIN, 11f));
            g2.drawString(n.label, n.x - 30, n.y + 26);
        }

        g2.dispose();
    }

    private static String pair(String a, String b) { return a + "→" + b; }

    interface PathSupplier { List<String> currentPath(); }
}
