import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Double;
import java.io.FileWriter;
import java.util.*;
import java.io.File;

/*
Remodel bias: Angular acceleration for bias at each time step is a
Gaussian centered at alpha times the value of the current bias. Find
alpha.
 */

public class findAlpha
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

    LinkedList <Pair< Double >> [] data;

    public static void main(String[] args)
    {
	findAlpha reader = new findAlpha();
	reader.getData();
    }

    public double gaussian(double x)
    {
	return Math.exp(-x*x / 2) / Math.sqrt(2*Math.PI);
    }

    public double gaussian(double x, double mu, double sigma)
    {
	return gaussian((x-mu)/sigma)/sigma;
    }

    public double getLL(int index, double param)
    {
	double ll = 0.0;

	Iterator<Pair<Double>> it = data[index].iterator();

	while(it.hasNext()){
	    Pair<Double> pair = it.next();

	    ll += Math.log(gaussian(pair.alpha,param*pair.omega,30.0));
	}

	return ll;
    }

    public void getData()
    {
	File pwd = new File(new File(".").getAbsolutePath());

	String files[] = new String[3];

	try{
	    files = new String[]{pwd.getCanonicalPath()+"/data/pimu_1",
				 pwd.getCanonicalPath()+"/data/pimu_2",
				 pwd.getCanonicalPath()+"/data/pimu_3"};
	}
	catch (IOException e) {
	    e.printStackTrace();
	}

	BufferedReader br = null;
	String line = "";
	String delim = ",";

	double[] val = new double[6]; //raw readings
	double[] dval = new double[6]; //differential of raw readings
	double[] w = new double[4]; //angular velocity
	double[] a = new double[4]; //angular acceleration

	for (String file: files){
	    data = new LinkedList[2];

	    for(int j=0; j<data.length; j++){
		data[j] = new LinkedList<Pair<Double>>();
	    }
	    
	    try{
		System.out.println("Parsing file:"+file);

		file = file + ".csv";

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

		    data[0].add(new Pair<Double>(w[0],a[0]));
		    data[1].add(new Pair<Double>(w[1],a[1]));

		    System.arraycopy(val, 0, prev, 0, val.length);
		    System.arraycopy(w, 0, prev, val.length, w.length);

		    count ++;

		    if (count%3.0E6 == 0)
			System.out.println("Line number:"+count);

		}

		System.out.println("Finished file, read "+count+" lines");

		//TODO: Insert call to function that finds alpha here
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
