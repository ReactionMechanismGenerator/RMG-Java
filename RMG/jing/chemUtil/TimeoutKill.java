/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jing.chemUtil;

/**
 *
 * gmagoon 12/4/08: cf.  Ch. 16 of "Ivor Horton's Beginning Java 2: JDK 5 Edition", esp. p. 727-728
 * this waits for the a specified amount of time, then kills the process (if it has not been interrupted) 
 */
public class TimeoutKill extends Thread{
    private Process process;
    private long delay;
    private String msg;
    
    public TimeoutKill(Process p_process, String p_msg, long p_delay){
        process=p_process;
        delay = p_delay;
        msg=p_msg;
        setDaemon(true);//this shouldn't really be relevant, but it is probably good practice to set it to true
    }
    
    public void run(){
        try {
            sleep(delay);
            process.destroy();
            //display message indicating the process has been killed
            System.out.println(msg);
        } catch (InterruptedException ex) {
            //do nothing since we expect to interrupt as typical behavior 
        }
        
        
    }
   

}
