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
    double[] time;

    int[] mean; //the estimate if we simply took a mean of the readings

    int N; //number of observations
    int Ns; //number of sensors
    int Nw; //number of possible 'omega' values

    int resW = 1; //resolution of state space
    
    double[] Tw; //transition matrix for omega
    double[][][] Tb; //transition matrix for each bias

    int b_lim = 40; //maximum allowed deviation of bias from the mean
    int db_lim = 25; //maximum change in bias allowed over a single time instant

    int[] b_means; //initial means for the biases

    int w_lim; //limits of the omega
    int[] states;

    double[][] V; //matrix to store log-likelihood
    int[][] B; //matrix to store most likely previous state

    //int[] seq; //best path sequence
    //int[] seq_backtrack; //best path by backtracking

    double[] w_model; //our model for the angular velocity
    double[][] b_model; //our model for the bias

    private Random rand;

    boolean DEBUG = false;

    double scale = 131.0;

    DecimalFormat df = new DecimalFormat("###.####");
    String prefix = "/home/dev/Documents/RAIG/data/";
    String file = prefix+"results";

    BufferedWriter logger = null;

    Viterbi(int N, int Ns, double[] w_model, double[][] b_model, long seed){

	if(seed == -1) 
	    rand = new Random();
	else 
	    rand = new Random(seed);

	this.N = N;
	
	this.Ns = Ns;

	this.w_model = w_model;
	this.b_model = b_model;

	this.obs = new int[this.Ns][this.N];
	this.bias = new int[this.Ns][this.N];
	this.time = new double[this.N];
	this.b_means = new int[this.Ns]; //the mean value of the bias for each sensor

	int buffer = 2000; // dirty hack

	readBias reader = new readBias(N+buffer, Ns);
	int[][] biases = reader.getBias();
	this.time = reader.getTime(); //the timestamps

	//use first 'buffer' values to get the mean, and the rest as data
	
	for(int s=0; s<this.Ns; s++){
	    for(int i=0; i<this.N+buffer; i++){
		if(i>=buffer)
		    this.bias[s][i-buffer] = biases[s][i];
		else
		    b_means[s] += biases[s][i];
	    }
	    b_means[s] /= buffer;
	}

	this.generateOmega(this.N); //generate omega values

	this.mean = new int[N]; //store the naive mean for baseline comparison

	if(bias[0].length != this.omega.length){
	    System.out.println("Number of observations don't match!!");
	    System.exit(0);
	}

	for(int sens=0; sens<this.Ns; sens++){
	    for(int t=0; t<this.N; t++){
		this.obs[sens][t] = this.bias[sens][t] + this.omega[t];
	    }
	}

	this.Nw = 2*(int)(Math.ceil(this.w_lim/this.resW)) + 1;
	this.Tw = new double[this.Nw];
	this.Tb = new double[this.Ns][this.b_lim*2+1][this.db_lim*2+1];

	this.states = new int[this.Nw];

	for(int i=0; i<this.Nw; i++)
	    this.states[i] = -this.w_lim + (i*this.resW);

	// probability of transitioning from 'r' to 'c'
	for(int r=0; r<this.Nw; r++)
	    this.Tw[r] = this.logProb((double)r*this.resW, this.w_model[0],
				      this.w_model[1]*scale);

	// cache transition probabilities for bias
	for(int s=0; s<this.Ns; s++){
	    for(int i=0; i<Tb[s].length; i++){
		int b = i-b_lim;
		double mu = b*b_model[s%3][0];
		double sig = b_model[s%3][1];
		for(int j=0; j<Tb[s][i].length; j++){
		    int db = j-db_lim;
		    this.Tb[s][i][j] = this.logProb((double)db, mu, sig);
		}
	    }
	}

	this.V = new double[2][this.Nw]; //store the log probabilities
	this.B = new int[this.N][this.Nw]; //store the best path for backtracking(can remove)

	for(int i=0; i<Nw; i++)
	    this.V[0][i] = this.logProb(this.states[i], this.omega[0], 0.5);
	//this.w_model[1]*scale);

	//this.seq = new int[this.N];
	//this.seq_backtrack = new int[this.N];

	System.out.println("Generated "+this.Nw+" states");
	System.out.println("Using "+this.Ns+" sensors for "+
			   this.N+" observations at a resolution of "+this.resW);

	System.out.println("bias= "+bias.length+" x "+bias[0].length);
	System.out.println("obs= "+obs.length+" x "+obs[0].length);
	System.out.println("omega= "+omega.length);

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
	//return logProb(fin-init, b_model[0], b_model[1]);
	return logProb(fin-init, 0.0, 7.0);
    }

    public double getBiasTrans(int sens, int init, int fin){
	
	int b_idx = (init-b_means[sens])-(-b_lim);
	int db_idx = (fin-init) - (-db_lim);

	//System.out.println("Indices for"+init+", "+fin+" are "+b_idx+", "+db_idx);
	
	if(b_idx<0 || b_idx>(2*b_lim) || db_idx<0 || db_idx>(2*db_lim))
	    return Double.NEGATIVE_INFINITY;
	else
	    return Tb[sens][b_idx][db_idx];
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

			p += getBiasTrans(sens, bi, bf);
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
		out = time[t]+","+df.format(omega[t])+","+df.format(seq[t])+
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

	if(args.length < 2){
	    System.out.println("Usage: [no. of observations] [no. of sensors] [seed]");
	    System.exit(0);
	}

	int N = Integer.parseInt(args[0]);
	int Ns = Integer.parseInt(args[1]);

	double[] w_model = new double[]{0.0, 0.1};
	double[][] b_model = new double[][]{{-0.9, 2.1}, {-0.91, 2.25}, {-0.95, 2.45}};
	//double[][] b_model = new double[][]{0.0, 7.0};

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
