package raig;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.File;
import java.io.IOException;
import java.lang.Double;
import java.util.ArrayList;

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
	String line = "";
	String delim = ", ";

	for(int i=0; i<files.length; i++){
	    try{
		String file = files[i];
		String file_output = file + "_mean.txt";
		file = file + ".txt";

		System.out.println("Processing file "+file);

		read = new BufferedReader(new FileReader(file));
		writer = new BufferedWriter(new FileWriter(file_output));

		String[] initial = (read.readLine()).split(delim);

		double prev_t = Double.parseDouble(initial[0])/1.0E6; 
		double t0 = prev_t;

		int count = 0;
		int numread = 0;
		int numwrote = 0;

		double sum_t = 0.0;
		double sum_r = 0.0;

		ArrayList<Double[]> data = new ArrayList<Double[]>();

		while((line = read.readLine()) != null){
		    numread++ ;
		    
		    String[] row = line.split(delim);

		    double t = Double.parseDouble(row[0])/1.0E6;
		    double dyaw = Double.parseDouble(row[1]);

		    double dt = t-prev_t;
		    dt += Math.abs(dt)<2.0 ? 0 : rollover;

		    count++;
		    t0 += dt;
		    
		    sum_t += t0;
		    sum_r += dyaw;

		    prev_t = t;

		    if(count == window){
			sum_r /= (count*scale);
			sum_t /= (count);

			data.add(new Double[]{sum_t, sum_r});

			String out = Double.toString(sum_t)+","+Double.toString(sum_r);
			writer.write(out);
			writer.newLine();

			numwrote++;
			
			count = 0;
			sum_r = sum_t = 0.0;
		    }
		}

		System.out.println("Read "+numread+" lines");
		System.out.println("Wrote "+numwrote+" lines");
		System.out.println("Arraylist size is "+data.size());
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
			writer.flush();
			writer.close();
		    } catch(IOException e) {
			e.printStackTrace();
		    }
		}
	    }
	}
    }

    public static void main(String[] args){
	
	String[] files = new String[]{"../data/bias0","../data/bias1", "../data/bias2"};

	parse r = new parse(files, 11);

	r.windowFitler();
    }
}
