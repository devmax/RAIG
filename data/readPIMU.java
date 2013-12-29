import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Double;
import java.io.FileWriter;

public class readPIMU
{
    public static void main(String[] args)
    {
	readPIMU reader = new readPIMU();
	reader.clean();
    }

    public void clean()
    {

	String files[] = {"/home/dev/Documents/RAIG/data/pimu_1",
			  "/home/dev/Documents/RAIG/data/pimu_2",
			  "/home/dev/Documents/RAIG/data/pimu_3"};

	BufferedReader br = null;
	FileWriter writer = null;
	String line = "";
	String delim = ",";

	double[] val = new double[9];

	for (String file: files){
	    try{
		System.out.println("Parsing file:"+file);

		String filtered = file + "_filtered.csv";
		file = file + ".csv";

		br = new BufferedReader(new FileReader(file));
		writer = new FileWriter(filtered);

		line = br.readLine(); //header
		String[] initial = (br.readLine()).split(delim); //First line
	    
		double[] init = new double[6];
		init[0] = Double.parseDouble(initial[0]); //time
		init[1] = Double.parseDouble(initial[1]); //pimu time
		init[2] = Double.parseDouble(initial[2]); // yaw 1
		init[3] = Double.parseDouble(initial[5]); // roll
		init[4] = Double.parseDouble(initial[6]); // yaw 2
		init[5] = Double.parseDouble(initial[9]); // pitch

		double ascale = 1.0/16384.0;

		int count = 0;
		writer.append("Time,Pimu Time,yaw1,roll,yaw2,pitch,ax,ay,az\n");

		while ((line = br.readLine()) != null){
		    String[] row = line.split(delim);

		    val[0] = (Double.parseDouble(row[0])-init[0])/1.0E6; //time
		    val[1] = (Double.parseDouble(row[1])-init[1])/1.0E6; //pimu time
		    val[2] = (Double.parseDouble(row[2])-init[2])/1.0E6 % 360; //yaw 1
		    val[3] = -(Double.parseDouble(row[5])-init[3])/1.0E6 % 360; //roll
		    val[4] = (Double.parseDouble(row[6])-init[4])/1.0E6 % 360; //yaw 2
		    val[5] = (Double.parseDouble(row[9])-init[5])/1.0E6 % 360; //pitch
		    val[6] = (Double.parseDouble(row[10]))*ascale; //ax
		    val[7] = (Double.parseDouble(row[11]))*ascale; //ay
		    val[8] = (Double.parseDouble(row[12]))*ascale; //az

		    for (int i=0; i<=7; i++){
			writer.append(Double.toString(val[i]));
			writer.append(',');
		    }
		    
		    writer.append(Double.toString(val[8]));
		    writer.append('\n');
		    count ++;

		    if (count%10000 == 0)
			System.out.println("Line number:"+count);
		}
		writer.flush();
		writer.close();

		System.out.println("Finished file");
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
