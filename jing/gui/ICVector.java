package jing.gui;

public class ICVector {
    public int Index = -1;
    public String Name;
    public String Concentration;
    public String Unit;
    public String React;
    public String AdjList;

    public ICVector(int ID, String name1, String Conc, String name2, String name3, String name4) {
    	this.Index = ID;
    	this.Name = name1;
    	this.Concentration = Conc;
    	this.Unit = name2;
    	this.React = name3;
    	this.AdjList = name4;
    }
}