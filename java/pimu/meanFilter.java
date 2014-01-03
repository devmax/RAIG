package pimu;

import java.util.LinkedList;
import java.lang.Double;
import java.util.ListIterator;
import java.util.Arrays;
import java.util.Scanner;

public class meanFilter
{
    int windowSize;
    LinkedList<Double> raw;
    LinkedList<Double> filtered;

    public meanFilter(int windowSize, LinkedList<Double> raw){
	this.windowSize = windowSize;
	this.raw = raw;

	filtered = new LinkedList<Double>();
    }

    public void filter(){
	int numel = 2*windowSize + 1;
	
	ListIterator<Double> left = raw.listIterator(0);
	ListIterator<Double> right = raw.listIterator(0);

	double sum=0;
	int index=0;
	while(index++<=windowSize){
	    sum +=right.next(); 
	}
	
	index=0;
	while(index<windowSize){
	    filtered.add(sum/(windowSize+index+1));
	    sum +=right.next();
	    index++;
	}

	while(right.hasNext()){
	    filtered.add(sum/numel);
	    sum +=right.next();
	    sum -=left.next();
	    index++;
	}
	
	index = raw.size()-index-1;
	while(index>=0){
	    filtered.add(sum/(windowSize+index+1));
	    sum -=left.next();
	    index--;
	}

	assert(raw.size() == filtered.size()) : "Sizes do not match!";
    }

    public LinkedList<Double> getFiltered(){
	filter();

	return filtered;
    }

    public static void main(String[] args){
	Scanner scan = new Scanner(System.in);

	int windowSize = scan.nextInt();
	LinkedList<Double> input = new LinkedList<Double>();
	LinkedList<Double> filtered = new LinkedList<Double>();
	double num;

	while((num=scan.nextDouble())!=-999.0){
	    input.add(num);
	}

	meanFilter mf = new meanFilter(windowSize,input);

	System.out.println("Input=\n"+Arrays.toString(input.toArray()));
	System.out.println("Result=\n"+Arrays.toString(mf.getFiltered().toArray()));
    }

}
