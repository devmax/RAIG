package raig;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.File;
import java.io.IOException;
import java.lang.Double;
import java.util.List;
import java.util.LinkedList;
import java.util.ListIterator;
import java.lang.Integer;

public class readBias
{
    String prefix = "/home/dev/Documents/RAIG/data/";

    /*
    String files[] = {"../../../../data/big0aa", "../../../../data/big1aa", "../../../../data/big2aa",
		      "../../../../data/big0ab", "../../../../data/big1ab", "../../../../data/big2ab",
		      "../../../../data/big0ac", "../../../../data/big1ac", "../../../../data/big2ac",
		      "../../../../data/big0ad", "../../../../data/big1ad", "../../../../data/big2ad",
		      "../../../../data/big0ae", "../../../../data/big1ae", "../../../../data/big2ae",
		      "../../../../data/big0af", "../../../../data/big1af", "../../../../data/big2af",
		      "../../../../data/big0ag", "../../../../data/big1ag", "../../../../data/big2ag" };
    */

    String[] files;

    int[][] bias;

    int N;
    int Ns;

    readBias(int N, int Ns){
	this.N = N;
	this.Ns = Ns;

	bias = new int[this.Ns][this.N];

	this.files = new String[this.Ns];

	String[] idx = {"a", "b", "c", "d", "e", "f", "g"};
	String name = "big";

	int lim = 3;

	System.out.print("Reading bias from: ");

	for(int i=0; i<Ns; i++){
	    files[i] = this.prefix+name+Integer.toString(i%lim)+idx[0]+idx[i/lim];
	    System.out.print(files[i]+", ");
	}

	System.out.println();
	
    }

    public int[][] getBias(){

	BufferedReader read = null;
	String line = "";
	String delim = ",";

	for(int i=0; i<Ns; i++){
	    try{
		String file = files[i];

		read = new BufferedReader(new FileReader(file));

		int numread = 0;

		while((line=read.readLine()) != null && numread<N){

		    String[] row = line.split(delim);
		    
		    bias[i][numread++] = (int)Double.parseDouble(row[1]);
		}

		if(numread != N)
		    System.out.println("Read only "+numread+" points");
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
		    } catch(IOException e) {
			e.printStackTrace();
		    }
		}
	    }
	}

	return bias;
    }
}
