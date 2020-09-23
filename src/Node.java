import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Node {

    int id;
    List<Link> links;

    public Node(int index){
        id = index;
        links = new ArrayList<>();
    }

    public int getId() {
        return id;
    }

    public void addLink(Link link) {
        links.add(link);
    }

    public boolean exceedLimit(int connectionLimit){
        return links.size() > connectionLimit;
    }

    public void removeRandomLink(Random rand) {
        int i = rand.nextInt(links.size());
        // Call method in link to remove itself from both nodes
        links.get(i).disconnect();
    }

    public void removeLink(Link link) {
        links.remove(link);
    }

    public Node[] checkNodesConnected(Node[] nodesConnected){
        // If this node is not already connected
        if(nodesConnected[id] != this){
            // Add this node to those connected
            nodesConnected[id] = this;
            // Recursively add the linked nodes
            for (Link l : links) {
                Node otherNode = l.getOtherNode(id);
                nodesConnected = otherNode.checkNodesConnected(nodesConnected);
            }
        }

        return nodesConnected;
    }

    public void addNewRandomLink(Node[] connectedNodes, Random rand, int linkLimit) {
        // Get a random node from those connected
        int i;
        do {
            i = rand.nextInt(connectedNodes.length);
        } while (connectedNodes[i] == null);
        // If the number of links on this node or the found node are at the limit
        if(this.exceedLimit(linkLimit - 1)){
            this.removeRandomLink(rand);
        }
        // Proceed to also check the found node
        if(connectedNodes[i].exceedLimit(linkLimit - 1)){
            connectedNodes[i].removeRandomLink(rand);
        }
        // Create a link between the found node and this node.
        new Link(connectedNodes[i], this);
    }
}
