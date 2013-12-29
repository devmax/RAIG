import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Double;
import java.io.FileWriter;

public class diff
{
    public static void main(String[] args)
    {
	diff reader = new diff();
	reader.clean();
    }

    public void clean()
    {
	String files[] = {"/Users/Dev/RAIG/data/pimu_1",
			  "/Users/Dev/RAIG/data/pimu_2",
			  "/Users/Dev/RAIG/data/pimu_3"};

	BufferedReader br = null;
	FileWriter writer = null;
	String line = "";
	String delim = ",";

	double[] val = new double[6]; //raw readings
	double[] dval = new double[6]; //differential of raw readings
	double[] w = new double[4]; //angular velocity
	double[] a = new double[4]; //angular acceleration

	for (String file: files){
	    try{
		System.out.println("Parsing file:"+file);

		String filtered = file + "_diff.csv";
		file = file + ".csv";

		br = new BufferedReader(new FileReader(file));
		writer = new FileWriter(filtered);

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
		writer.append("time,yaw1,yaw2,roll,pitch,d_yaw1,d_yaw2,d_roll,d_pitch,a_yaw1,a_yaw2,a_roll,a_pitch\n");

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

		    System.arraycopy(val, 0, prev, 0, val.length);
		    System.arraycopy(w, 0, prev, val.length, w.length);

		    writer.append(Double.toString(val[0]));
		    writer.append(',');

		    for (int i=2; i<=5; i++){
			writer.append(Double.toString(val[i]));
			writer.append(',');
		    }

		    for (int i=0; i<=3; i++){
			writer.append(Double.toString(w[i]));
			writer.append(',');
		    }

		    for (int i=0; i<=3; i++){
			writer.append(Double.toString(a[i]));
			if(i!=3)
			    writer.append(',');
		    }

		    writer.append('\n');
		    
		    count ++;


		    if (count%3.0E6 == 0)
			System.out.println("Line number:"+count);

		}
		writer.flush();
		writer.close();

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
