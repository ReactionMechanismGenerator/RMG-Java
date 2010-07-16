package jing.chem;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

//6/3/09 gmagoon: molFile class to store and perform operations on molFiles
//may eventually want to have a super class threeDGeom and say "public class molFile extends threeDGeom {"
public class molFile {
    //attributes
    private String name;
    private String directory;
    
    //constructors
    //Constructor 0: trivial, empty constructor
    public molFile(){
    }
    //Constructor 1: construct a molFile object while writing a file with text determined by chemGraph
    public molFile(String p_name, String p_dir, ChemGraph p_chemGraph){
        //assign attributes
        name=p_name;
        directory=p_dir;
        //first, create a string
        String molFileString = p_chemGraph.generateMolFileString();//generate the string
        //write the string to file (code taken from Mike's code in Species.java
        File mf = null;
        try {
            mf = new File(directory+"/"+name+".mol");
            FileWriter fw = new FileWriter(mf);
            fw.write(molFileString);
            fw.close();
        } catch (IOException e) {
            String err = "Error writing " + directory+name+".mol";
            err += e.toString();
            System.out.println(err);
        }
        
    }
    //Constructor 2: construct a molFile object that already exists; the object will be a pointer to the location of the file; it is assumed that the molFile already exists at the path specified
    public molFile(String p_name, String p_dir){
        //assign attributes
        name=p_name;
        directory=p_dir;
    }

    public String getPath(){
        return directory+"/"+name+".mol";
    }

    public String getCrudePath(){
        return directory+"/"+name+".cmol";
    }
    
    public String getName(){
        return name;
    }
    
    public String getDirectory(){
        return directory;
    }
    
    //eventually, we can have operations for returning a molefile string by reading in a moleFile
    
}
