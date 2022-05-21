import it.unisa.dia.gas.jpbc.Element;

import java.util.Arrays;

public class Node {
    //Gate is represented by two numbers (T, n), n is the number of child nodes, and t is the threshold value
    //Null if it is a leaf node
    public int[] gate;
    //Children refers to child nodes. This field is the index list of child nodes
    //Null if it is a leaf node
    public int[] children;
    //Att represents the attribute value,
    //If it is an internal node, this field is null
    public int att;
    // the secret of the node
    public Element secretShare;

    // to decide whether  the node can be recovered
    public boolean valid;
    //way to construct an internal node
    public Node(int[] gate, int[] children){
        this.gate = gate;
        this.children = children;
    }
    // way to construct a leaf node
    public Node(int att){
        this.att = att;
    }
    public boolean isLeaf() {
        return this.children == null;
    }
    @Override
    public String toString() {
        if (this.isLeaf()){
            return Integer.toString(this.att);
        }
        else {
            return Arrays.toString(this.gate);
        }
    }
}
