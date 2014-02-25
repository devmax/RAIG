package raig;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.lang.ArrayIndexOutOfBoundsException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.File;
import java.io.IOException;
import java.lang.Double;
import java.util.Random;
import java.util.Arrays;
import java.text.DecimalFormat;

class Viterbi
{
    int[] omega;
    int[][] bias;
    int[][] obs;

    int[] mean; //the estimate if we simply took a mean of the readings

    int N; //number of observations
    int Ns; //number of sensors
    int Nw; //number of possible 'omega' values

    int resW = 1; //resolution of state space
    
    double[] Tw; //transition matrix

    int w_lim;
    int[] states;

    double[][] V; //matrix to store log-likelihood
    int[][] B; //matrix to store most likely previous state

    int[] seq; //best path sequence
    int[] seq_backtrack; //best path by backtracking

    double[] v_error; //error between truth and Viterbi estimate
    double[] m_error; //error between truth and mean estimate

    double[] w_model; //our model for the angular velocity
    double[] b_model; //our model for the bias

    private Random rand;

    boolean DEBUG = false;

    double scale = 131.0;

    Viterbi(int N, int Ns, double[] w_model, double[] b_model, long seed){

	if(seed == -1) 
	    rand = new Random();
	else 
	    rand = new Random(seed);

	for(int i=0; i<1150; i++)
	    rand.nextGaussian();

	this.N = N;
	
	this.Ns = Ns;

	this.w_model = w_model;
	this.b_model = b_model;

	this.obs = new int[this.Ns][this.N];

	readBias reader = new readBias(N, Ns);

	this.bias = reader.getBias();

	this.generateOmega(this.N);

	this.mean = new int[N];

	if(bias[0].length != this.omega.length){
	    System.out.println("Number of observations don't match");
	    System.exit(0);
	}

	for(int sens=0; sens<this.Ns; sens++){
	    for(int t=0; t<this.N; t++){
		this.obs[sens][t] = this.bias[sens][t] + this.omega[t];
	    }
	}

	this.Nw = 2*(int)(Math.ceil(this.w_lim/this.resW)) + 1;
	System.out.println("Generated "+this.Nw+" states");
	System.out.println("Using "+this.Ns+" sensors for "+
			   this.N+" observations at a resolution of "+this.resW);

	System.out.println("bias= "+bias.length+" x "+bias[0].length);
	System.out.println("obs= "+obs.length+" x "+obs[0].length);
	System.out.println("omega= "+omega.length);

	this.Tw = new double[this.Nw];

	this.states = new int[this.Nw];

	for(int i=0; i<this.Nw; i++)
	    this.states[i] = -this.w_lim + (i*this.resW);

	// probability of transitioning from 'r' to 'c'
	for(int r=0; r<this.Nw; r++)
	    this.Tw[r] = this.logProb((double)r*this.resW, this.w_model[0],
				      Math.round(this.w_model[1]*scale));

	this.V = new double[2][this.Nw];
	this.B = new int[this.N][this.Nw];

	for(int i=0; i<Nw; i++)
	    this.V[0][i] = this.logProb(this.states[i], this.omega[0], 0.5);
	//this.w_model[1]*scale);

	this.seq = new int[this.N];
	this.seq_backtrack = new int[this.N];

	this.v_error = new double[this.N];
	this.m_error = new double[this.N];

	if(DEBUG){

	    System.out.println("Tw is:");
	    dispArray(Tw);

	    System.out.println("States are:");
	    dispArray(states);

	    System.out.println("Omega is:");
	    dispArray(omega);

	    System.out.println("Bias is:");
	    dispArray(bias);

	    System.out.println("Obs is:");
	    dispArray(obs);
	}

    }

    public static void dispArray(double[][] arr){
	for(int i=0; i<arr.length; i++)
	    dispArray(arr[i]);
    }

    public static void dispArray(double[] arr){
	System.out.print("[");

	for(int i=0; i<arr.length; i++)
	    System.out.format("%.3f, ", arr[i]);

	System.out.println("]");
    }

    public static void dispArray(int[][] arr){
	for(int i=0; i<arr.length; i++)
	    dispArray(arr[i]);
    }

    public static void dispArray(int[] arr){
	System.out.print("[");

	for(int i=0; i<arr.length; i++)
	    System.out.format("%d, ", arr[i]);

	System.out.println("]");
    }

    public int[] generateOmega(int N){
	omega = new int[N];
	w_lim = (int)Double.NEGATIVE_INFINITY;

	for(int t=1; t<N; t++){
	    omega[t] = omega[t-1] + (int)Math.round(getGaussian(this.w_model[0],
								this.w_model[1])*scale);
	    if(Math.abs(omega[t])>w_lim){
		w_lim = Math.abs(omega[t]);
            }
	}

        w_lim = Math.abs(w_lim)+10;

	System.out.println("Omega limit is: "+w_lim);

	return omega;
    }

    private double getGaussian(double mu, double var){
	return mu + rand.nextGaussian()*var;
    }

    public double logProb(double x, double mu, double sigma){
	return -(Math.log(sigma) + 0.5*Math.pow((x-mu)/sigma, 2));
    }

    public double getBiasTrans(int init, int fin){
	return logProb(fin-init, b_model[0], b_model[1]);
    }

    public void findSequence(){

	double p_ml;
	int ml = 0;

	double p_max;
	int s_max;

	double p;
        int t;

	for(t=1; t<N; t++){ //loop over time
	    p_ml = Double.NEGATIVE_INFINITY;
	    ml = -1;

	    if(t%100 == 0)
		System.out.print("Time = "+t+"\n");

	    for(int i=0; i<Nw; i++){ //trying new states
		p_max = Double.NEGATIVE_INFINITY;
		s_max = -1;

		for(int j=0; j<Nw; j++){ //trying old states
		    p = V[0][j] + Tw[Math.abs(j-i)]; //initialize log-probability
		    for(int sens=0; sens<Ns; sens++){
			int bi = obs[sens][t-1]-states[j];
			int bf = obs[sens][t]-states[i];

			p += getBiasTrans(bi, bf);
		    }

		    if(p > p_max){
			p_max = p;
			s_max = j;
		    }
		}
		V[1][i] = p_max;
		B[t][i] = s_max;

		if(p_max > p_ml){
		    p_ml = p_max;
		    ml = i;
		}
	    }

	    seq[t] = states[ml];

	    if(DEBUG){
		System.out.println("V= ");
		dispArray(V);
	    }

	    System.arraycopy(V[1], 0, V[0], 0, V[1].length);

	    for(int sens=0; sens<Ns; sens++)
		mean[t] += obs[sens][t]-bias[sens][0];    

	    mean[t] /= Ns;

	    v_error[t] = (seq[t]-omega[t]);///scale;
	    m_error[t] = (mean[t]-omega[t]);///scale;

	    if(DEBUG){
		System.out.format("GT = %d ,", omega[t]);
		System.out.format("dw = %d ,", omega[t]-omega[t-1]);
		System.out.format("v_err = %.4f, ", v_error[t]);
		System.out.format("mu_err = %.4f\n", m_error[t]);
	    }
	}

	for(t=N-1; t>=0; t--){
	    seq_backtrack[t] = states[ml];
	    ml = B[t][ml];
	}
    }

    private void Log(){

	DecimalFormat df = new DecimalFormat("###.####");

        String prefix = "/home/dev/Documents/RAIG/data/";
	String file = prefix+"results";
	//String err = "../data/error";

	String out;

	BufferedWriter logger = null;
	//BufferedWriter error = null;

	try{
	    logger = new BufferedWriter(new FileWriter(file));

	    for(int t=0; t<N; t++){
		out = df.format(omega[t])+","+df.format(seq[t])+
		    ","+df.format(mean[t])+","+df.format(seq_backtrack[t]);
		logger.write(out);
		logger.newLine();

		/*
		out = df.format(v_error[t])+","+df.format(m_error[t]);
		error.write(out);
		error.newLine();
		*/

	    }
	}
	catch (IOException e) {
	    e.printStackTrace();
	}
	finally{
	    if(logger != null){
		try{
		    logger.flush();
		    logger.close();

		    /*
		    error.flush();
		    error.close();
		    */
		}catch (IOException e){
		    e.printStackTrace();
		}
	    }
	}
    }

    public static void main(String[] args){

	int N = Integer.parseInt(args[0]);
	int Ns = Integer.parseInt(args[1]);

	double[] w_model = new double[]{0.0, 0.1};
	double[] b_model = new double[]{0.0, 7.0};

	Viterbi v;

	if(args.length == 3){
	    System.out.println("Using seed "+args[2]);
	    v = new Viterbi(N, Ns, w_model, b_model, Long.parseLong(args[2]));
	}
	else
	    v = new Viterbi(N, Ns, w_model, b_model, -1);
	
	v.findSequence();
	v.Log();
    }
}
