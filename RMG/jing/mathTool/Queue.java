////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
//	RMG Team (rmg_dev@mit.edu)
//
//	Permission is hereby granted, free of charge, to any person obtaining a
//	copy of this software and associated documentation files (the "Software"),
//	to deal in the Software without restriction, including without limitation
//	the rights to use, copy, modify, merge, publish, distribute, sublicense,
//	and/or sell copies of the Software, and to permit persons to whom the
//	Software is furnished to do so, subject to the following conditions:
//
//	The above copyright notice and this permission notice shall be included in
//	all copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//	DEALINGS IN THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////



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

