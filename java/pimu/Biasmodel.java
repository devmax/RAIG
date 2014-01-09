package pimu;

import java.lang.Double;
import java.lang.ArithmeticException;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Arrays;
import java.lang.Integer;

import pimu.pimureader;

public class Biasmodel
{
    LinkedList<Double> omega;
    LinkedList<Double> alpha;

    double theta;
    pimureader reader;

    int skip;

    public Biasmodel(int limit,int skip, int index, String file)
    {
	reader=new pimureader(limit,skip,new int[]{index},new String[]{file});
	reader.filter(10);

	this.omega=reader.getOmega()[0][0];
	this.alpha=reader.getAlpha()[0][0];

	/*
	ListIterator<Double> w = this.omega.listIterator(150);
	ListIterator<Double> a = this.alpha.listIterator(150);

	int count=0;
	while(count++<15)
	    System.out.println("w:"+w.next()+", a:"+a.next());
	*/

    }

    public void computeTheta()
    {
	double num=0,den=0;
	ListIterator<Double> w=omega.listIterator(skip);
	ListIterator<Double> a=alpha.listIterator(skip);
	while(w.hasNext() && a.hasNext()){
	    double rate=w.next();
	    double accel=a.next();
	    num += rate*accel;
	    den += rate*rate;
	}
	theta=(num/den);
    }

    public double getTheta(){
	return theta;
    }

    public static void main(String[] args)
    {
	int limit=(int)11.0E6;
	int skip=25000;
	int step=(int)4.0E6;

	if(args.length>=1){
	    limit=Integer.parseInt(args[0]);
	}
	if(args.length>=2){
	    skip=Integer.parseInt(args[1]);
	}
	if(args.length>=3){
	    step=Integer.parseInt(args[2]);
	}

	int[] indices=new int[]{0,1};
	String[] files=new String[]{"pimu_1","pimu_2","pimu_3"};
	String[] axes=new String[]{"yaw_1","yaw_2","roll","pitch"};

	System.out.println();
	for(int j=0;j<files.length;j++){

	    System.out.println("Processing file "+files[j]);
	    for(int k=0;k<indices.length;k++){
		System.out.println("Processing axis "+axes[indices[k]]);

		LinkedList<Double> theta=new LinkedList<Double>();
		LinkedList<Integer> sizes=new LinkedList<Integer>();
		double total=0;
		System.out.println("*****************************************");

		for(int i=skip; i<limit; i+=step){
		    System.out.println("Processing lines "+i
				       +" to "+Math.min(i+step,limit));
		    Biasmodel armodel = new Biasmodel(Math.min(i+step,limit),
						      i,indices[k],files[j]);
		    armodel.computeTheta();
		    System.out.println(armodel.getTheta());
		    System.out.println("-------------------------------");

		    theta.add(armodel.getTheta());
		    sizes.add(armodel.omega.size());
		    total += armodel.omega.size();
		}

		double weighted=0;

		ListIterator<Double> param=theta.listIterator(0);
		ListIterator<Integer> s=sizes.listIterator(0);
		
		while(param.hasNext()){
		    weighted += (double)(s.next()/total)*param.next();
		}
		System.out.println("\nWEIGHTED estimate: "+weighted);
		System.out.println("*****************************************");
	    }
	}

	
    }
}
