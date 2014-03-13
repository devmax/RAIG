package raig;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.Double;
import java.util.Random;
import java.util.Arrays;
import java.text.DecimalFormat;

class Viterbi
{

    int w;
    int[] b;
    int[][] obs;
    double[][] time;
    int mean;
    double theta;

    int t;
    double dt;

    int N; //number of observations
    int Ns; //number of sensors
    int Nw; //number of possible 'omega' values

    int resW = 1; //resolution of state space
    
    double[] Tw; //transition matrix for omega
    double[][][] Tb; //transition matrix for each bias

    int b_lim = 40; //maximum allowed deviation of bias from the mean
    int db_lim = 25; //maximum change in bias allowed over a single time instant

    int[] b_means; //initial means for the biases

    int w_lim=650; //limits of the omega
    int[] states;

    double[][] V; //matrix to store log-likelihood
    double[][] yaw; //matrix to store all possible values of yaw so far

    double[] w_model; //our model for the angular velocity
    double[][] b_model; //our model for the bias

    private Random rand;

    double scale = 131.0;

    DecimalFormat df = new DecimalFormat("###.####");
    String prefix = "/home/dev/Documents/RAIG/data/";
    String delim=",";

    String[] files;
    String outfile = prefix+"results";

    BufferedWriter logger;
    BufferedReader reader[];

    Viterbi(String[] files, int N, double[] w_model, double[][] b_model, long seed){

	if(seed == -1) 
	    rand = new Random();
	else 
	    rand = new Random(seed);

	this.N = N;
	
	this.Ns = files.length;

	this.w_model = w_model;
	this.b_model = b_model;

	System.out.println("Reading from files: ");

	reader = new BufferedReader[this.Ns];

	try{
	    for(int i=0; i<Ns; i++){
		System.out.print(files[i]+" ");
		reader[i] = new BufferedReader(new FileReader(prefix+files[i]));
	    }
	}
	catch(FileNotFoundException e){
	    e.printStackTrace();
	}

	System.out.println();

	try{
	    logger = new BufferedWriter(new FileWriter(outfile));
	}
	catch (IOException e) {
	    e.printStackTrace();
	}

	this.b = new int[this.Ns];
	this.obs = new int[2][this.Ns];
	this.time = new double[2][this.Ns];

	this.Nw = 2*(int)(Math.ceil(this.w_lim/this.resW)) + 1;
	this.Tw = new double[this.Nw];
	this.Tb = new double[this.Ns][this.b_lim*2+1][this.db_lim*2+1];

	this.states = new int[this.Nw];
	this.V = new double[2][this.Nw]; //store the log probabilities
	this.yaw = new double[2][this.Nw];

	this.t = 0;
	this.dt = 0.0;
	
	this.b_means = new int[this.Ns]; //the mean value of the bias for each sensor

	int buffer = 2000; // dirty hack

	for(int i=0; i<buffer; i++){
	    this.nextBias();
	    for(int s=0; s<this.Ns; s++)
		b_means[s] += this.b[s];
	}

	for(int s=0; s<this.Ns; s++)
	    b_means[s] = Math.round(b_means[s]/buffer);

	System.out.println("Bias means are:"+Arrays.toString(b_means));

	for(int i=0; i<this.Nw; i++)
	    this.states[i] = -this.w_lim + (i*this.resW);

	// probability of transitioning from 'r' to 'c'
	for(int r=0; r<this.Nw; r++){
	    this.yaw[0][r] = this.yaw[1][r] = 0.0;
	    this.Tw[r] = this.logProb((double)r*this.resW, this.w_model[0],
				      this.w_model[1]*scale);
	}

	// cache transition probabilities for bias
	for(int s=0; s<this.Ns; s++){
	    for(int i=0; i<Tb[s].length; i++){
		int bias = i-b_lim;
		double mu = bias*b_model[s%3][0];
		double sig = b_model[s%3][1];
		for(int j=0; j<Tb[s][i].length; j++){
		    int db = j-db_lim;
		    this.Tb[s][i][j] = this.logProb((double)db, mu, sig);
		}
	    }
	}

	this.nextObs();
	System.arraycopy(obs[1], 0, obs[0], 0, obs[1].length);

	for(int i=0; i<Nw; i++)
	    this.V[0][i] = this.logProb(this.states[i], this.w, 0.5);
	//this.w_model[1]*scale);

	System.out.println("Generated "+this.Nw+" states");
	System.out.println("Using "+this.Ns+" sensors for "+
			   this.N+" observations at a resolution of "+this.resW);

    }

    public void nextObs(){

	int w_p = w;

	System.arraycopy(obs[1], 0, obs[0], 0, obs[1].length);
	System.arraycopy(time[1], 0, time[0], 0, time[1].length);

	nextBias();
	nextOmega();

	dt = 0.0;

	for(int i=0; i<Ns; i++)
	    obs[1][i] = b[i]+w;

	for(int i=0; i<Ns; i++){
	    if(time[1][i]-time[0][i] < 0){
		System.out.println("Bad timestamp detected! "+
				   time[0][i]+"->"+time[1][i]+
				   ",diff="+(time[1][i]-time[0][i]));
		break;
	    }
	    dt += time[1][i]-time[0][i];
	}

	dt /= Ns;

	theta += w_p*dt + 0.5*(w-w_p)*dt*dt;
    }

    public void nextBias(){
	String line;
	String[] row;

	try{
	    for(int i=0; i<Ns; i++){
		if((line = reader[i].readLine()) != null){
		    row = line.split(delim);
		    b[i] = (int)Double.parseDouble(row[1]);
		    time[1][i] = Double.parseDouble(row[0]);
		}
		else{
		    System.out.println("Nothing to read...");
		    N = 0;
		}
	    }
	}
	catch(IOException e){
	    e.printStackTrace();
	}
    }

    public void nextOmega(){

	if(t==0){
	    w=0;
	}
	else{

	    int dw = (int)Math.round(getGaussian(this.w_model[0], this.w_model[1])*scale);

	    w = Math.min(w_lim, Math.max(-w_lim, w+dw));
	}
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

	int estimate;

	String out;
	int p_mean;
	double mu_theta = 0.0;

	while(t<N){ //loop over time

	    nextObs();

	    p_ml = Double.NEGATIVE_INFINITY;
	    ml = -1;

	    if(t%10000 == 0)
		System.out.print("Time = "+t+"\n");

	    for(int i=0; i<Nw; i++){ //trying new states
		p_max = Double.NEGATIVE_INFINITY;
		s_max = -1;

		for(int j=0; j<Nw; j++){ //trying old states
		    p = V[0][j] + Tw[Math.abs(j-i)]; //initialize log-probability
		    for(int sens=0; sens<Ns; sens++){
			int bi = obs[0][sens]-states[j];
			int bf = obs[1][sens]-states[i];

			p += getBiasTrans(sens, bi, bf);
		    }

		    if(p > p_max){
			p_max = p;
			s_max = j;
		    }
		}
		V[1][i] = p_max;

		if(s_max != -1)
		    yaw[1][i] = yaw[0][s_max] + 
			dt*states[s_max] + 0.5*(states[i]-states[s_max])*dt*dt;

		if(p_max > p_ml){
		    p_ml = p_max;
		    ml = i;
		}
	    }

	    estimate = states[ml];
	    
	    System.arraycopy(V[1], 0, V[0], 0, V[1].length);
	    System.arraycopy(yaw[1], 0, yaw[0], 0, yaw[1].length);
	    p_mean = mean;
	    mean = 0;

	    for(int sens=0; sens<Ns; sens++)
		mean += obs[1][sens]-(b_means[sens]);

	    mean = Math.round((float)mean/Ns);
	     
	    mu_theta += dt*p_mean + 0.5*(mean-p_mean)*dt*dt;

	    out = Double.toString(dt)+","+
		Double.toString(w)+","+Double.toString(estimate)+","+Double.toString(mean)+","+
		Double.toString(theta/scale)+","+Double.toString(yaw[1][ml]/scale)+","+
		Double.toString(mu_theta/scale);

	    try{
		logger.write(out);
		logger.newLine();
	    }
	    catch(IOException e){
		e.printStackTrace();
	    }

	    ++t;
	}

	for(int i=0; i<Ns; i++){
	    if(reader[i] != null){
		try{
		    reader[i].close();
		}
		catch(IOException e){
		    e.printStackTrace();
		}
	    }
	}

	try{
	    logger.flush();
	    logger.close();
	}
	catch(IOException e){
	    e.printStackTrace();
	}

    }

    public static void main(String[] args){

	if(args.length < 2){
	    System.out.println("Usage: [no. of observations] [seed(-1 for none)] [file 1]...[file Ns]");
	    System.exit(0);
	}

	int N = Integer.parseInt(args[0]);

	long seed = Long.parseLong(args[1]);

	double[] w_model = new double[]{0.0, 0.1};
	double[][] b_model = new double[][]{{-0.9, 2.1}, {-0.91, 2.25}, {-0.95, 2.45}};
	//double[][] b_model = new double[][]{0.0, 7.0};

	String[] files = new String[args.length-2];

	for(int i=0; i<files.length; i++)
	    files[i] = args[i+2];

	Viterbi v = new Viterbi(files, N, w_model, b_model, seed);
	
	v.findSequence();
    }
}
