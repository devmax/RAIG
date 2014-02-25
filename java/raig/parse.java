package raig;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.File;
import java.io.IOException;
import java.lang.Double;
import java.util.List;
import java.util.LinkedList;
import java.util.ListIterator;

public class parse
{
    String files[];
    
    int window;
    public final static double scale = 131.0;

    public final static double rollover = Math.pow(2, 32)/1.0E6;

    public parse(String[] f, int w){
	this.files = f;
	this.window = w;
    }

    public void windowFitler(){
	
	BufferedReader read = null;
	BufferedWriter writer = null;
	BufferedWriter writer_reg = null;
	String line = "";
	String delim = ", ";

	boolean writeMean = false;

	for(int i=0; i<files.length; i++){
	    try{
		String file = files[i];

		if(writeMean){
		    String file_output = file + "_sub.txt"; // file to store the moving window averages
		    writer = new BufferedWriter(new FileWriter(file_output));
		}

		String file_reg = file + "_reg.txt"; // file to store the data for the OLS regression
		file = file + ".txt";

		System.out.println("Processing file "+file);

		read = new BufferedReader(new FileReader(file));
		writer_reg = new BufferedWriter(new FileWriter(file_reg));

		String[] initial = (read.readLine()).split(delim);

		double prev_t = Double.parseDouble(initial[0])/1.0E6; 
		double t0 = prev_t;

		int count = 0;
		int numread = 0;
		int numwrote = 0;

		double sum_t = 0.0;
		double sum_r = 0.0;

		List<Double[]> data = new LinkedList<Double[]>();

		while((line = read.readLine()) != null){
		    numread++ ;
		    
		    String[] row = line.split(delim);

		    double t = Double.parseDouble(row[0])/1.0E6;
		    double dyaw = Double.parseDouble(row[1]);

		    double dt = t-prev_t;
		    dt += Math.abs(dt)<2.0 ? 0 : rollover;

		    count++;

		    // this is the true new time, so just add dt to the previous time
		    t0 += dt;
		    
		    sum_t += t0;
		    sum_r += dyaw;

		    // remember the time we read at this instant for the future diff subtractions
		    prev_t = t;

		    if(count == window){
			sum_r /= (count*scale);
			sum_t /= (count);

			data.add(new Double[]{sum_t, sum_r});

			if(writeMean){
			    String out = Double.toString(sum_t)+","+Double.toString(sum_r);
			    writer.write(out);
			    writer.newLine();
			}

			numwrote++;
			
			count = 0;
			sum_r = sum_t = 0.0;
		    }
		}

		System.out.println("Read "+numread+" lines");
		System.out.println(numwrote+" lines after moving window average");


		numwrote = 0;

		ListIterator<Double[]> s_prev = data.listIterator(0);
		ListIterator<Double[]> s_next = data.listIterator(1);

		while(s_next.hasNext()){
		    Double[] s0 = s_prev.next();
		    Double[] s1 = s_next.next();

		    String out = Double.toString(s1[0]-s0[0])+","+Double.toString(s0[1])+
			","+Double.toString(s1[1]-s0[1]); // dt, b_0, b_1-b_0
		    writer_reg.write(out);
		    writer_reg.newLine();

		    numwrote++;
		}

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
			if(writeMean){
			    writer.flush();
			    writer.close();
			}
			writer_reg.flush();
			writer_reg.close();
		    } catch(IOException e) {
			e.printStackTrace();
		    }
		}
	    }
	}
    }

    public static void main(String[] args){
	
	String[] files = new String[]{"../data/big0","../data/big1", "../data/big2"};

	parse r = new parse(files, 11);

	r.windowFitler();
    }
}
