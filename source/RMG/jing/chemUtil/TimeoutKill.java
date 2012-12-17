// //////////////////////////////////////////////////////////////////////////////
//
// RMG - Reaction Mechanism Generator
//
// Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
// RMG Team (rmg_dev@mit.edu)
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// //////////////////////////////////////////////////////////////////////////////
package jing.chemUtil;

/**
 * gmagoon 12/4/08: cf. Ch. 16 of "Ivor Horton's Beginning Java 2: JDK 5 Edition", esp. p. 727-728 this waits for the a
 * specified amount of time, then kills the process (if it has not been interrupted)
 */
public class TimeoutKill extends Thread {
    private Process process;
    private long delay;
    private String msg;

    public TimeoutKill(Process p_process, String p_msg, long p_delay) {
        process = p_process;
        delay = p_delay;
        msg = p_msg;
        setDaemon(true);// this shouldn't really be relevant, but it is probably good practice to set it to true
    }

    public void run() {
        try {
            sleep(delay);
            process.destroy();
            // display message indicating the process has been killed
            System.out.println(msg);
        } catch (InterruptedException ex) {
            // do nothing since we expect to interrupt as typical behavior
        }
    }
}
