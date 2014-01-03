package pimu;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Double;
import java.lang.ArithmeticException;
import java.io.FileWriter;
import java.util.*;
import java.io.File;

import org.math.plot.*;

import pimu.readPIMU;

/*
Remodel bias: Angular acceleration for bias at each time step is a
Gaussian centered at alpha times the value of the current bias. Find
alpha.
 */

public class findAlpha
{
    LinkedList<Double>[][] omega;
    LinkedList<Double>[][] alpha;
    double initial;

    readPIMU reader;

    double learning_rate;
    double eg;
    double ex;
    int maxiter;

    int[] indices;

    public findAlpha(){
	indices = new int[]{0};
	
	reader = new readPIMU(1000, indices, new String[]{"pimu_1","pimu_2","pimu_3"});
	reader.filter(5);

	omega=reader.getOmega();
	alpha=reader.getAlpha();

	initial = 0.1;

	eg = 0.1E-3;
	ex = 0.1E-3;
	maxiter = 20;

	learning_rate = 0.1;
    }

    public void setInitAlpha(String arg)
    {
	this.initial = Double.parseDouble(arg);
    }

    public static void main(String[] args)
    {
	findAlpha descent = new findAlpha();

	if(args.length == 1)
	    descent.setInitAlpha(args[0]);
	else
	    System.out.println("Using default initial value of "+descent.initial+" for alpha");

	/*
	for(int i=0; i<descent.omega.length; i++){
	    System.out.println("\nPerforming gradient descent on file "+(i+1));
	    for(int j=0; j<descent.indices.length; j++){
		System.out.println("Index="+(j+1));
		descent.gradientDescent(i, j);
	    }
	}
	*/

    }

    public double gaussian(double x)
    {
	return Math.exp(-x*x / 2) / Math.sqrt(2*Math.PI);
    }

    public double gaussian(double x, double mu, double sigma)
    {
	return gaussian((x-mu)/sigma)/sigma;
    }

    public double getLL(int file, int index, double param)
    {
	double ll = 0.0;
	double srate = reader.getSRate(file);

	ListIterator<Double> w=omega[file][index].listIterator();
	ListIterator<Double> a=alpha[file][index].listIterator();

	int count = 0;
	while(w.hasNext()){
	    count ++;
	    
	    double rate = w.next();
	    double accel = a.next();
	    double mu = -(param/srate)*rate;

	    try{

		double pdf = Math.log(gaussian(accel,mu,30.0));
		if(pdf == Double.NEGATIVE_INFINITY)
		    throw new ArithmeticException("INFINITY EXCEPTION: x="+
						  accel+",mu="+mu+",w="+rate);
		ll += pdf;
	    }	    
	    catch(ArithmeticException e){
		System.out.println(e.getMessage());
	    }
		
	}
	return ll;
    }

    public void gradientDescent(int file, int idx)
    {
	double lrate = learning_rate;

	boolean converged = false;
	double param = initial;

	double prevll = getLL(file, idx, param);

	System.out.println("\nIteration 1"+", For alpha="+param+":");
	System.out.println("LL="+prevll);

	double newparam, ll, der;

	int count = 1;

	while (!converged){
	    count ++;

	    newparam = (1-lrate)*param;
	    
	    ll = getLL(file, idx, param);

	    der = (ll-prevll);//(newparam-param);

	    System.out.println("\nIteration "+count+", For alpha="+newparam+":");
	    System.out.println("LL="+ll+", Change="+der);

	    if(count > maxiter){
		System.out.println("Exceeded maximum number of iterations, stopping descent");
		break;
	    }
	    param = newparam;
	    prevll = ll;
	}
    }


}
