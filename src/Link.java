public class Link {
    Node _n1;
    Node _n2;

    public Link(Node n1, Node n2){
        _n1 = n1;
        _n2 = n2;

        // Add this Link object to each node
        _n1.addLink(this);
        _n2.addLink(this);
    }

    public void disconnect() {
        _n1.removeLink(this);
        _n2.removeLink(this);
    }

    public Node getOtherNode(int currentNodeID){
        if(_n1.getId() == currentNodeID){
            return _n2;
        }
        else{
            return _n1;
        }
    }
}
