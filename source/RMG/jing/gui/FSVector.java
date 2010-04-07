package jing.gui;

public class FSVector {
    public int Index = -1;
    public String Name;
    public String AdjList;

    public FSVector(int ID, String name, String graph) {
    	this.Index = ID;
    	this.Name = name;
    	this.AdjList = graph;
    }
}