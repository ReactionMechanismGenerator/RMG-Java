package jing.gui;

public class SMVector {
    public int Index = -1;
    public String Name;
    public String React;
    public String Location;

    public SMVector(int ID, String name1, String name2, String name3) {
    	this.Index = ID;
    	this.Name = name1;
    	this.React = name2;
    	this.Location = name3;
    }
}