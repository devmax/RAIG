package raig;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.lang.NumberFormatException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.File;
import java.io.IOException;
import java.lang.Double;
import java.lang.Integer;
import java.util.List;
import java.util.LinkedList;
import java.util.ListIterator;

/*
  Compute the histogram of the joint distribution of (bias, dbias)
*/
public class parse
{
    String files[];
    
    int window;
    public final static double scale = 1;//131.0;

    int maxCount;

    int bRange;
    int bMin;

    int dbMin;

    int[][] hist;

    //public final static double rollover = Math.pow(2, 32)/1.0E6;

    public parse(String[] f, int w, int num){
	this.files = f;
	this.window = w;

	if (num == -1)
	    this.maxCount = Integer.MAX_VALUE; //this is a hack, and
					       //will fail when
					       //num. of samples >
					       //2^31
	else
	    this.maxCount = num;

	bRange = 35;
	dbMin = -25;
    }

    public void windowFitler(){
	
	BufferedReader read = null;
	BufferedWriter writer_reg = null;
	BufferedWriter writer_hist = null;
	String line = "";
	String delim = ", ";

	boolean writeReg = true;

	for(int i=0; i<files.length; i++){

	    hist = new int[2*bRange + 1][-2*dbMin + 1];

	    try{
		String file = files[i];

		if(writeReg){
		    String file_reg = file + "_reg.txt"; // file to store the data for the OLS regression
		    writer_reg = new BufferedWriter(new FileWriter(file_reg));		    
		}

		String file_hist = file + "_hist.txt"; // file to store the frequency of each (dbias,bias) pair
		writer_hist = new BufferedWriter(new FileWriter(file_hist));

		file = file + ".txt";

		read = new BufferedReader(new FileReader(file));

		System.out.println("Processing file "+file);

		int numread = 0;
		int numwrote = 0;

		double sum_t = 0.0;
		double sum_r = 0.0;

		double min_r = Double.POSITIVE_INFINITY;
		double max_r = Double.NEGATIVE_INFINITY;

		double dbias; 
		double[] prev = new double[3]; // prev mean, min, max
		boolean inited = false;

		while((line = read.readLine()) != null && numread<maxCount){
		    numread++ ;
		    
		    String[] row = line.split(delim);

		    double t = Double.parseDouble(row[0])/1.0E6;
		    double dyaw = Double.parseDouble(row[1]);
		    
		    sum_t += t;
		    sum_r += dyaw;

		    min_r = min_r > dyaw? dyaw: min_r; // minimum in window
		    max_r = max_r < dyaw? dyaw: max_r; // maximum in window

		    if(numread%window == 0){

			sum_r = Math.round(sum_r/(window*scale));
			sum_t /= window;

			if(!inited){
			    bMin = (int)sum_r-bRange;

			    inited = true;
			}
			else{
			    dbias = sum_r - prev[0];

			    try{
				hist[(int)prev[0]-bMin][(int)dbias-dbMin]++; // increase histogram count				
			    }
			    catch(ArrayIndexOutOfBoundsException e){
				System.out.println("Trying to access "+prev[0]+","+dbias);
			    }


			    if(writeReg){
				String out = sum_t+","+prev[0]+","+dbias+","+prev[1]+","+prev[2]; 
				// t, mean_b, d(mean_b), min_b, max_b

				writer_reg.write(out);
				writer_reg.newLine();
			    }

			    numwrote++;
			}
			
			prev[0] = sum_r;
			prev[1] = min_r;
			prev[2] = max_r;
			sum_r = sum_t = 0.0;
			min_r = Double.POSITIVE_INFINITY;
			max_r = Double.NEGATIVE_INFINITY;
		    }
		}

		String lims = "Bias:["+bMin+","+(bMin+2*bRange)+"]\n"+
		    "Dbias:["+dbMin+","+(-dbMin)+"]";

		int row,col;

		writer_hist.write(lims);
		writer_hist.newLine();

		System.out.println("Histogram limits are:\n"+lims);

		for(row=0; row<hist.length; row++){
		    for(col=0; col<hist[row].length-1; col++){
			writer_hist.write(hist[row][col]+",");
		    }
		    writer_hist.write(hist[row][col]+"\n");
		}

		System.out.println("Read "+numread+" lines");
		System.out.println(numwrote+" lines after moving window average");
		System.out.println("Wrote "+numwrote+" lines for regression");
	    }
	    catch(FileNotFoundException e) {
		e.printStackTrace();
	    }
	    catch(IOException e) {
		e.printStackTrace();
	    }
	    finally{
		if (read != null) {
		    try{
			read.close();

			if(writeReg){
			    writer_reg.flush();
			    writer_reg.close();
			}

			writer_hist.flush();
			writer_hist.close();

		    } catch(IOException e) {
			e.printStackTrace();
		    }
		}
	    }
	}
    }

    public static void main(String[] args){

	String prefix = "/home/dev/Documents/RAIG/data/";
	String[] files = new String[args.length-2];

	int window=11, lim=1000;

	try{
	    window = Integer.parseInt(args[0]);
	    lim = Integer.parseInt(args[1]);
	}
	catch(NumberFormatException e){
	    System.out.println("Usage: [window size] [number of samples(-1 for all)] "+
			       "[file 1] .... [file N]");

	    System.exit(0);
	}


	for(int i=2; i<args.length; i++){
	    files[i-2] = prefix+args[i];
	}

	System.out.println("Window size = "+window);
	parse r = new parse(files, window, lim);

	r.windowFitler();
    }
}
