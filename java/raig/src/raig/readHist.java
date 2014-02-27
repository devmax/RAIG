package raig;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Integer;
import java.lang.Double;
import java.util.List;
import java.util.Arrays;

public class readHist
{
    String file;

    double[][] hist;

    readHist(String file){
	this.file = file;

	hist = new double[71][51];
    }

    public double[][] getHist(){

	return hist;
    }

    public double[] getParsed(String[] row){
	double[] count = new double[row.length];

	double sum = 0;
	for(int i=0; i<row.length; i++){
	    count[i] = Double.parseDouble(row[i]);
	    sum += count[i];
	}

	if(sum>0)
	    for(int i=0; i<count.length; i++)
		count[i] /= sum;

	return count;
    }

    public void read(){
	
	BufferedReader read = null;

	String line="";
	String delim = ",";

	try{
	    read = new BufferedReader(new FileReader(file));

	    read.readLine(); //bias range
	    read.readLine(); //dbias range

	    int r = 0;

	    while((line = read.readLine()) != null){
		
		hist[r++] = getParsed(line.split(delim));
	    }

	    System.out.println("Read "+hist.length+"x"+hist[0].length+" matrix");
	}
	catch(FileNotFoundException e){
	    e.printStackTrace();
	}
	catch(IOException e) {
	    e.printStackTrace();
	}
	finally{
	    if(read!=null){
		try{
		    read.close();
		}
		catch(IOException e){
		    e.printStackTrace();
		}
	    }
	}
    }

    public static void main(String[] args){
	
	String file = "/home/dev/Documents/RAIG/data/train1_hist.txt";

	readHist reader = new readHist(file);

	reader.read();

	double[][] h = reader.getHist();

	for(int i=0; i<h.length; i++){
	    for(int j=0; j<h[i].length; j++)
		System.out.print(h[i][j]+",");
	    System.out.println("\b ");
	}
    }
}
