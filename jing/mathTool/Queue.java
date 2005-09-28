//!********************************************************************************
//!
//!    RMG: Reaction Mechanism Generator                                            
//!
//!    Copyright: Jing Song, MIT, 2002, all rights reserved
//!     
//!    Author's Contact: jingsong@mit.edu
//!
//!    Restrictions:
//!    (1) RMG is only for non-commercial distribution; commercial usage
//!        must require other written permission.
//!    (2) Redistributions of RMG must retain the above copyright
//!        notice, this list of conditions and the following disclaimer.
//!    (3) The end-user documentation included with the redistribution,
//!        if any, must include the following acknowledgment:
//!        "This product includes software RMG developed by Jing Song, MIT."
//!        Alternately, this acknowledgment may appear in the software itself,
//!        if and wherever such third-party acknowledgments normally appear.
//!  
//!    RMG IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED 
//!    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
//!    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
//!    DISCLAIMED.  IN NO EVENT SHALL JING SONG BE LIABLE FOR  
//!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
//!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
//!    OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;  
//!    OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  
//!    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
//!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
//!    THE USE OF RMG, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//! 
//!******************************************************************************



package jing.mathTool;


import java.util.*;

//## package jing::mathTool 

//----------------------------------------------------------------------------
// jing\mathTool\Queue.java                                                                  
//----------------------------------------------------------------------------

//## class Queue 
public class Queue {
    
    private int first = -1;		//## attribute first 
    
    private int last = -1;		//## attribute last 
    
    private int size = 0;		//## attribute size 
    
    private Object [] storage;		//## attribute storage 
    
    
    // Constructors
    
    //## operation Queue(int,Collection) 
    public  Queue(int p_maxSize, Collection p_collection) {
        //#[ operation Queue(int,Collection) 
        int presentSize = p_collection.size();
        storage = new Object[Math.max(presentSize*2,p_maxSize)];
        
        Iterator iter = p_collection.iterator();
        while (iter.hasNext()) {
        	Object o = (Object)iter.next();
        	enqueue(o);
        }
        
        
        
        //#]
    }
    //## operation Queue(int) 
    public  Queue(int p_maxSize) {
        //#[ operation Queue(int) 
        if (p_maxSize < 0) p_maxSize = 20;
        storage = new Object[p_maxSize];
        
        
        //#]
    }
    public  Queue() {
    }
    
    //## operation clear() 
    public void clear() {
        //#[ operation clear() 
        first = last = -1;
        size = 0;
        
        
        //#]
    }
    
    //## operation dequeue() 
    public Object dequeue() {
        //#[ operation dequeue() 
        if (isEmpty()) throw new EmptyQueueException();
        
        Object o = storage[first];
         
        if (first == last) last = first = -1;
        else if (first == size -1) first = 0;
        else first++;
        
        return o;
        	
        
        	  
        
        
        //#]
    }
    
    //## operation enqueue(Object) 
    public void enqueue(Object p_object) {
        //#[ operation enqueue(Object) 
        if (isFull()) throw new FullQueueException();
        
        if (last == size -1 || last == -1) {
        	last = 0;
        	storage[last] = p_object;
        	if (first == -1) first = 0;
        }
        else{
        	storage[++last]= p_object;
        }
        return;
        	
        
        
        
        //#]
    }
    
    //## operation isEmpty() 
    public boolean isEmpty() {
        //#[ operation isEmpty() 
        return (first == -1);
        //#]
    }
    
    //## operation isFull() 
    public boolean isFull() {
        //#[ operation isFull() 
        return (first == 0 && last == size -1) || (first == last + 1);
        
        
        //#]
    }
    
    //## operation peek() 
    public Object peek() {
        //#[ operation peek() 
        if (first == -1) throw new EmptyQueueException();
        return storage[first];
        //#]
    }
    
}
/*********************************************************************
	File Path	: RMG\RMG\jing\mathTool\Queue.java
*********************************************************************/

