package pimu;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Double;
import java.util.*;

public class readPIMU
{
    public class Pair<A> 
    {
	private final A omega;
	private final A alpha;

	public Pair(A omega, A alpha){
	    this.omega = omega;
	    this.alpha = alpha;
	}

	public A getOmega(){
	    return omega;
	}

	public A getAlpha(){
	    return alpha;
	}

    }

    String files[]; 
    boolean read;

    int limit;
    int[] indices; //0:yaw1, 1:yaw2, 2:roll, 3:pitch

    LinkedList <Pair< Double >> [][] data;
    double srate[];

    public readPIMU(int limit, int[] indices, String[] files)
    {
	this.limit = limit;
	this.indices = indices;

	this.files = new String[files.length];

	for(int i=0; i<files.length; i++){
	    this.files[i] = "/home/dev/Documents/RAIG/data/"+files[i]+".csv";
	}

	this.data = new LinkedList[this.files.length][this.indices.length];

	for(int i=0; i<this.files.length; i++){
	    for(int j=0; j<this.indices.length; j++){
		this.data[i][j] = new LinkedList<Pair<Double>>();
	    }
	}

	this.read = false;
	this.srate = new double[this.files.length];
    }

    public readPIMU(int[] indices, String[] files)
    {
	this(-1, indices, files);
    }

    public LinkedList<Pair<Double>>[][] getData()
    {
	if(!read)
	    readData();

	return data;
    }

    public double getSRate(int file)
    {
	return srate[file];
    }

    public boolean checkRollover(double prev, double current)
    {
	if((prev+current) < 10.0 && (prev+current) > -10.0 && Math.abs(prev-current) > 500)
	    return true;
	else
	    return false;
    }

    public void readData()
    {
	read = true;

	BufferedReader br = null;
	String line = "";
	String delim = ",";


	double[] val = new double[6]; //raw readings
	double[] dval = new double[6]; //differential of raw readings
	double[] w = new double[4]; //angular velocity
	double[] a = new double[4]; //angular acceleration

	for (int i=0; i<files.length; i++){
	    
	    try{
		String file = files[i];
		System.out.println("\nParsing file:"+file);

		br = new BufferedReader(new FileReader(file));

		line = br.readLine(); //header
		String[] initial = (br.readLine()).split(delim); //First line
	    
		double[] prev = new double[10];
		prev[0] = Double.parseDouble(initial[0]); //time
		prev[1] = Double.parseDouble(initial[1]); //pimu time
		prev[2] = Double.parseDouble(initial[2]); // yaw 1
		prev[3] = -Double.parseDouble(initial[5]); // roll
		prev[4] = Double.parseDouble(initial[6]); // yaw 2
		prev[5] = -Double.parseDouble(initial[9]); // pitch

		int count = 0;

		while ((line = br.readLine()) != null){
		    String[] row = line.split(delim);

		    val[0] = Double.parseDouble(row[0]); //time
		    val[1] = Double.parseDouble(row[1]); //pimu time
		    val[2] = Double.parseDouble(row[2]); //yaw 1
		    val[3] = -Double.parseDouble(row[5]); //roll
		    val[4] = Double.parseDouble(row[6]); //yaw 2
		    val[5] = -Double.parseDouble(row[9]); //pitch

		    if(checkRollover(prev[2],val[2])){
			System.out.println("Rollover detected in yaw1:"+
					   prev[2]+"->"+val[2]);
			val[2] += prev[2];
		    }
		    if(checkRollover(prev[3],val[3])){
			System.out.println("Rollover detected in roll:"+
					   prev[3]+"->"+val[3]);
			val[3] += prev[3];
		    }
		    if(checkRollover(prev[4],val[4])){
			System.out.println("Rollover detected in yaw2:"+
					   prev[4]+"->"+val[4]);
			val[4] += prev[4];
		    }
		    if(checkRollover(prev[5],val[5])){
			System.out.println("Rollover detected in pitch:"+
					   prev[5]+"->"+val[5]);
			val[5] += prev[5];
		    }

		    dval[0] = (val[0]-prev[0])/1.0E6; //time
		    dval[1] = (val[1]-prev[1])/1.0E6; //pimu time
		    dval[2] = (val[2]-prev[2])/1.0E6; //yaw 1
		    dval[3] = (val[3]-prev[3])/1.0E6; //roll
		    dval[4] = (val[4]-prev[4])/1.0E6; //yaw 2
		    dval[5] = (val[5]-prev[5])/1.0E6; //pitch

		    w[0] = dval[2]/dval[0]; //dyaw 1
		    w[1] = dval[4]/dval[0]; //dyaw 2
		    w[2] = dval[3]/dval[0]; //droll
		    w[3] = dval[5]/dval[0]; //dpitch

		    a[0] = (w[0]-prev[6])/dval[0]; //d2yaw 1
		    a[1] = (w[1]-prev[7])/dval[0]; //d2yaw 2
		    a[2] = (w[2]-prev[8])/dval[0]; //d2roll
		    a[3] = (w[3]-prev[9])/dval[0]; //d2pitch
		    
		    for(int j=0; j<indices.length; j++){
			data[i][j].add(new Pair<Double>(w[indices[j]],a[indices[j]]));
		    }

		    System.arraycopy(val, 0, prev, 0, val.length);
		    System.arraycopy(w, 0, prev, val.length, w.length);

		    count ++;

		    if (count%3.0E6 == 0)
			System.out.println("Line number:"+count);

		    if (count == limit){
			//System.out.println("Reached limit, stop reading");
			break;
		    }

		}
		srate[i] = dval[0];
		System.out.println("Finished file, read "+count+" lines");

	    }
	    catch (FileNotFoundException e) {
		e.printStackTrace();
	    }
	    catch (IOException e) {
		e.printStackTrace();
	    }
	    finally {
		if (br != null) {
		    try {
			br.close();
		    } catch (IOException e) {
			e.printStackTrace();
		    }
		}
	    }
	}
    }
}
