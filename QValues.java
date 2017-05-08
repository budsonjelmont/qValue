
import java.io.*;

/**
 *	Program QValue:
 *		This program takes an XML document, output by Filemaker, gets all of the experiment data
 *		from the document and reformats it, and uses this reformatted data for input to a 
 *		Java-to-R interface. R calculates the p values and then q values for this data, and the 
 *		list of q values it outputs is put into a format that Filemaker can import to mark the 
 *		significance of each experiment contained therein.
 *
 *	Class QValues:
 *		This class does all the work described above except the XML parsing and the Java-to-R interfacing. 
 *		For the former it instantiates the XMLToDataArray2 class, and for the latter it instantiates the
 *		TalkToR class.
 * 
 */

public class QValues {

	int num_timepoints = -1;			// These are set to their default values, so that if the command line is
	int num_replicates = -1;			// looked through and no argument is there to set a certain variable, that
	int lf_weight = 1;					// variable will keep its default value. (Default values of -1 mean that
	int silac_weight_1 = -1;			// by default this variable will not be considered; in the cases of
										// "timepoints" and "replicates," if these still = -1 even after
	int silac_weight_2 = -1;			// reading the command line, there will be an error.)
	
	int min_replicates = 2;
	int ref_tpt = 0;		
	int _unpaired=-1;
	int BHfdr=0;
	String path_to_file, file_name;
	double min_peak_area = 0;
	String path_to_wait;
	boolean log_lf = false;
	boolean log_silac = false;
	boolean _error = false;
	private TalkToR r_talker;
   
	/**
	 * @param possible command line arguments:
	 * 
	 *  (args[0] does not get a variable name, but it should be the file path+name that Filemaker output)
	 *  path_to_wait = the file path+name of the text file that gives a boolean whether the program is still running
	 * 	num_timepoints = the (max) number of timepoints at which replicates were tested
	 * 	num_replicates = the (max) number of replicates tested
	 * 	lf_weight = the SILAC weight selected for label free comparison
	 * 	silac_weight_1 = the SILAC weight in the numerator of the SILAC ratio (SILAC ratio.above)
	 * 	silac_weight_2 = the SILAC weight in the denominator of the SILAC ratio (SILAC ratio.below)
	 * 	min_replicates = the minimum number of replicates desired with which to perform the significance tests
	 *  ref_tpt = the label-free timepoint to which to compare the other timepoints: can be either "min" or "max"
	 *  min_peak_area = the minimum peak-area value that should be considered (indicates all values below are "noise")
	 *  log_lf = whether or not to take the log (base 2) of the label free peakareas before sending to R
	 *  log_silac = whether or not to take the log (base 2) of the silac peakareas before sending to R
	 *  _error = whether or not an error file should be printed if something goes wrong (for testing purposes)
	 *  bh = whether to control FDR using the Benjamini-Hochberg method instead of Storey
	 * All possible command line arguments:
	 * 
	 * 	filename pathToWait timepoints=x replicates=x label_free=x silac_numerator=x silac_denominator=x min_replicates=x min_peak_area=x bh=x min/max logLF logSILAC error
	 * 
	 * -> filename, pathToWait, timepoints, replicates datatype are NECESSARY, but these and the other paremeters, except filename
	 *    and pathToWait, can be put on the command line in any order
	 * 
	 */

	public static void main(String[] args) {
		if(args.length<4 || args.length>14){
			System.out.println("Error: Invalid number of parameters1");
			return;
		} 
		QValues q=new QValues(); 
		
		String[] forPTF = args[0].split("\\\\");
 		int endLen = forPTF[forPTF.length-1].length();
 		q.path_to_file = args[0].substring(0, args[0].length()-endLen);		 		 		
 		String name=forPTF[forPTF.length-1];
  		String[] name2=name.split("\\.");
  		q.file_name=q.path_to_file+name2[0];
  		
		q.path_to_wait = args[1];		
		q.wait(true);	
		// set the program as busy = true as it calculates the q values
		// set the instance variables to their non-default values, if they were input
		for(int i = 2; i < args.length; i++){		// start at the 2nd argument -> this is where the "=x" args start
			String[] list = args[i].split("=");						// look through each argument for what it's setting and the intended value
			if(list.length>1){										// if a parameter is "param=", with no input number, don't
				if(list[0].equals("timepoints")){					// set that parameter
					q.num_timepoints = Integer.parseInt(list[1]);
				}
				if(list[0].equals("replicates")){
					q.num_replicates = Integer.parseInt(list[1]);
				}
								
				if(list[0].equals("unpaired")){
					q._unpaired = Integer.parseInt(list[1]);
				}
				if(list[0].equals("label_free")){
					q.lf_weight = Integer.parseInt(list[1]);
				}
				if(list[0].equals("silac_numerator")){
					q.silac_weight_1 = Integer.parseInt(list[1]);
				}
				if(list[0].equals("silac_denominator")){
					q.silac_weight_2 = Integer.parseInt(list[1]);
				}
				if(list[0].equals("min_replicates")){
					q.min_replicates = Integer.parseInt(list[1]);
				}
				if(list[0].equals("min_peak_area")){
					q.min_peak_area = Double.parseDouble(list[1]);
				}
				if(list[0].equals("bh")){
					q.BHfdr = Integer.parseInt(list[1]);
				}
			}// if the parameter doesn't get input with an equals sign, will be checked for here
			if(list[0].equals("error")){							// if the parameter is "error," set error to true
				q._error = true;
			}
			if(list[0].equals("max")){								// by default this is set = 0, which is the "min" setting
				q.ref_tpt = 1;										// so for input of "max," change val to 1
			}
			if(list[0].equals("logLF")){
				q.log_lf = true;
			}
			if(list[0].equals("logSILAC")){
				q.log_silac = true;
			}
		}
		// check that timepoints and replicates were set-- these are NECESSARY to run the program
		if(q.num_timepoints == -1 || q.num_replicates == -1){
			System.out.println("Error: Invalid parameters2");
			return;
		}
		// if you have timepoints, replicates (and assuming the filename was first) then continue!	
		boolean hasSilac = false;
		if(q.silac_weight_1!=-1 && q.silac_weight_2!=-1){	// if the silac weights are (both) given
			hasSilac = true;		// (if only one given, still can't calculate q values for silac)
		}
		System.out.println("Gathering data...");		
		
		// parse the XML file from Filemaker and ignore any values below minimum peak area input by user
		XMLToDataArray2 x = new XMLToDataArray2(q.num_timepoints, q.min_peak_area);
		double[][] data = x.parseFilemakerData(args[0]);
//		q.printOutputArray(data); //print parsed XML from FM for debugging
		
		try {
			q.r_talker = new TalkToR();
			q.r_talker.openR();	
			System.out.println("Performing statistical calculations on label-free data...");
			// use the parsed XML data to get the label-free q-values
			double[][] labelFree = q.makeLabelFreeList(data);		
			double[][] SILAC = new double[0][0];
			
			if(hasSilac){	// if the silac weights are given, use them and get SILAC q-values
				System.out.println("Performing statistical calculations on SILAC data...");
				SILAC = q.makeSILACList(data);
			}
			
			// now labelFree and SILAC are the 2D arrays q-values for these, formatted per tpt per experiment
			System.out.println("Reformatting Q Values for FileMaker...");
			double[][] listForFilemaker = q.consolidate(labelFree, SILAC, data, hasSilac);
			System.out.println("Printing Filemaker input...");
			
			q.printFilemakerInput(args[0], listForFilemaker);
			System.out.println("Done!");
			q.wait(false);									// set as busy = false now that the program is done
		} catch (BadDataException e) {
			q.reportError(e.getMessage());
		}
	}

	/**
	 * makeLabelFreeList:
	 *   This method looks through the array of data taken from the XML document, finds which
	 *   timepoint is the reference, then gets the sums of the data and the standard deviations
	 *   in order to compare each non-reference timepoint to the reference. If "logLF" was entered
	 *   on the command line, then this method gets the sums and standard deviations of the logs
	 *   (base 2) of the data. Then it arranges them in the proper order for R: a 6-column table
	 *   with the sums of the each experiment's data, the number of data points summed, (presumably
	 *   these to get averages) and their standard deviations. If this method finds an index in the
	 *   data array with no data (denoted by Double.NaN) it does not include that in the sum
	 *   or standard deviation calculation. It also does not include a data point in the sum or
	 *   standard deviation if it comes from calcuating the sum, average, or standard deviation of
	 *   data without the minimum number of data points specified by the user (or the default, 
	 *   which is 2).
	 * 	
	 * @param array, the array of data obtained by parsing the xml document
	 * @return output, a table arranged in the proper way for Filemaker to import q-values
	 * 		   based on experiment row and timepoint column
	 */
	
	double[][] makeLabelFreeList(double[][] array) throws BadDataException{
		if(array == null){
			throw new BadDataException("No data with which to make the R input", path_to_wait, _error);
		}
		System.out.println("\t" + "(Organizing the data)");
		// array for R is a 2D array: one dimension will be of fixed length (the double[] to be entered per
		// timepoint). The other is an ArrayList to accommodate the uncertainty of how many timepoints will
		// have statistics calculated: at most (tpts - 1) per experiment (since 1 tpt must be the reference)
		// but some may not be calculated due to lack of data, so can't know beforehand
		java.util.ArrayList<double[]> arrayForR = new java.util.ArrayList<double[]>();
		// output will be the q-values for label-free: num_tpts cols (ref tpt will be blank) by num experiments rows
		double[][] output = new double[num_timepoints][array[0].length];
		//System.out.println(array.length);
		for(int i = 0; i < array[0].length; i++){				// for each row in the excel data array
			boolean hasData = false;
			for(int j = 1; j < array.length; j++){				// j = 1 -- > start after the counter
				if(!(Double.isNaN(array[j][i]))){
					hasData = true;								// has data if even one index has a value
				}
			}
			if(hasData){
				double[] sums = new double[num_timepoints];
				double[] numSummed = new double[num_timepoints];
				// Get the sums of each timepoint's data
				for(int j = 0; j < sums.length; j++){
					for(int k = 0; k < num_replicates; k++){	// num_replicates = highest possible number of replicates	
						int lookingAt = (k*num_timepoints) + (j+1) + ((num_timepoints*num_replicates)*(lf_weight-1));
						//the horiz. ind//this is which	//this is how	// this adds the number of results that come before
						//you're looking//timept u want	//offset the	// the results you want depending on which SILAC 
						//@ in the table//in the table	//starting pt is// you want get significance values for
						if(!(Double.isNaN(array[lookingAt][i]))){
							if(log_lf){									// if command line directed to take the log
								sums[j] = sums[j] + (Math.log(array[lookingAt][i])/Math.log(2));
								// add log of peakarea to the sums: Math.log is ln (as in, base e) but lnx/ln2 = log (base 2) x
							}
							else{
								sums[j] = sums[j] + array[lookingAt][i];
							}
							numSummed[j]++;
						}		
					}	// if it's not an actual value, sums[j] = numSummed[j] = 0
				}
				// Find the minimum or maximum average and which timepoint it is
				int refIndex = 0;
				if(ref_tpt==0){								// if it's set to compare to the min peak area tpt
					double minAvg = Double.POSITIVE_INFINITY;
					for(int j = 0; j < sums.length; j++){
						if(numSummed[j] > 1){
							if((sums[j]/numSummed[j]) < minAvg){
								minAvg = (sums[j]/numSummed[j]);
								refIndex = j;
							}
							// -if there isn't a sum, don't bother looking because it's not the minimum
							// -if there's only one number summed, it won't have a standard deviation so
							//  its significance can't be ascertained (min_replicates taken into account later)
						}					
					}
				}
				else if(ref_tpt==1){						// if it's set to compare to the max peak area tpt
					double maxAvg = Double.NEGATIVE_INFINITY;
					for(int j = 0; j < sums.length; j++){
						if(numSummed[j] > 1){
							if((sums[j]/numSummed[j]) > maxAvg){
								maxAvg = (sums[j]/numSummed[j]);
								refIndex = j;
							}
						}					
					}
				}
				double[] stDevs = this.calcStDevs(array, i, sums, numSummed, lf_weight);
				// now have the reference tpt, the sums of the peakareas, the number of replicates per tpt
				//  --> can now set up to calculate p-values (and later q values)!
				for(int j = 0; j < sums.length; j++){
					if(j != refIndex && numSummed[j]>=min_replicates){
						double[] thisTpt = new double[6];
						// if it's not the reference (if reference, can't compare it to itself)
						// or the timepoint didn't have enough data to get an accurate standard deviation
						thisTpt[0] = sums[refIndex];
						thisTpt[1] = sums[j];
						thisTpt[2] = numSummed[refIndex];
						thisTpt[3] = numSummed[j];
						thisTpt[4] = stDevs[refIndex];
						thisTpt[5] = stDevs[j];	
						arrayForR.add(thisTpt);
					}
					else{
						// if it is the reference, or it didn't have enough data
						output[j][i] = -1;
					}
				} 
			}
			else{	// if this row doesn't have data, make num_tpts rows (for each timepoint) and put -1 in each spot
				for(int j = 0; j < num_timepoints; j++){
					output[j][i] = -1;
				}

			}
		}
		// transfer the data from the ArrayList format to the double array format (for TalktoR's "get q values" method)
		double[][] tableForR = new double[6][arrayForR.size()];
		for(int j = 0; j < arrayForR.size(); j++){
			double[] curr = arrayForR.get(j);
			for(int i = 0; i < 6; i++){
				tableForR[i][j] = curr[i];
			}
		}
				
		System.out.println("\t" + "(Getting small-sample t-test p-values)");
		String file1=file_name+"_labelFree_pHist.png";
		String file2=file_name+"_labelFree_qplot.png";
		double[] pVals = r_talker.getSmallSamplePVals(tableForR);
		
		////Code block below prints computed p values
//		java.lang.StringBuffer buffer = new java.lang.StringBuffer();
//		// Build the table
//	   for(int col = 0; col < pVals.length; col++){
//			buffer.append(pVals[col] + "\r\n");					
//		}
//		String pValOutput = buffer.toString();
//		try{
//			String[] forPTF = path_to_wait.split("\\\\");
//			int endLen = forPTF[forPTF.length-1].length();
//			String pathToFile = path_to_wait.substring(0, path_to_wait.length()-endLen) + "_pVals for troubleshooting.txt";
//			// Create file 
//			FileWriter fstream = new FileWriter(pathToFile);
//			BufferedWriter out = new BufferedWriter(fstream);
//			out.write(pValOutput);
//			//Close the output stream
//			out.close();
//		}catch (Exception e){//Catch exception if any
//			System.err.println("Error: " + e.getMessage());
//		}
		////End pVal printing
		
		// pVals is a double array: need in ArrayList format for R Talker's "get qvalues" method
		java.util.ArrayList<Double> pValsForR = new java.util.ArrayList<Double>();
		for(int i = 0; i < pVals.length; i++){
			pValsForR.add(pVals[i]);
		}
	
		// now in correct format, so send it to R to get the q-values!
		System.out.println("\t" + "(Correcting p-values for multiple tests: getting q-values)");
		double[] qVals = r_talker.getQValues(pValsForR,file1,file2,BHfdr);
		// now, need to return q-values in the format they'll go into Filemaker
		// - all non-calculable q-value spots (non-calculable timepoints) already have placeholders in output
		//	 (-1s), and there are as many blank (not -1) spaces in 2D output as there are qvalues in qVals
		//	 list : so put the q-values from R one by one into the blank spaces -> because of the way the
		//	 p-values were sent to R (while keeping track of timepoint spaces in output), this results in 
		//	 q-values inserted into the right places in output
		int qVRow = 0;
		for(int j = 0; j < array[0].length; j++){
			for(int i = 0; i < num_timepoints; i++){
				if(output[i][j] > -1){				// if the row in table for R denoting this timepoint was calculable
					output[i][j] = qVals[qVRow];		// then put the next available qvalue into this timepoint in output
					qVRow++;							// get ready to put in next available qvalue from qval list
				}										// (works because qvals only made for calculable data in tableForR)
				// no else case, because if timepoint wasn't calculable, there's already a placeholder (-1) in that spot
			}
		}
		
		return output;
	}
	
	/**
	 * silacCalculations:
	 *    This method looks through the array of data taken from the XML document, and for each 
	 *   timepoint inserts into a table its replicate measurements of a chosen SILAC weight and
	 *   the corresponding replicate measurements of a different chosen SILAC weight. The table,
	 *   a 2D double array, looks like this:		
	 *    This table is then sent to R and			[rep 1 weight 1 for tpt 0] [rep 1 weight 2 for tpt 0]
	 *   reproduced, after which the paired		  	[rep 2 weight 1 for tpt 0] [rep 2 weight 2 for tpt 0]
	 *   t-test is performed on it and the			[rep 3 weight 1 for tpt 0] [rep 3 weight 2 for tpt 0]
	 *   p-value for likelihood of a difference		[rep 4 weight 1 for tpt 0] [rep 4 weight 2 for tpt 0]
	 *   between the SILAC weights returned. This	[rep 5 weight 1 for tpt 0] [rep 5 weight 2 for tpt 0]
	 *   is done for each timepoint in each row of	
	 *   the table. The t-test p-values are returned here, and then sent again to TalkToR to calculate
	 *   the q-values for the SILAC comparisons.
	 *    If this method finds an index in the data array with no data (denoted by Double.NaN) it
	 *   inserts a -1 into the table; TalkToR handles this missing data and all possible consequences
	 *   therein when reproducing the table in R or receiving the result of a statistical call to R. 
	 * 	  If SILAC data is provided in the XML document by the user, this method will be called and
	 *   will return the SILAC "half" of the total results table, a 2D String array of each timepoint's
	 * 	 q-value.
	 * 
	 * @param array, the array of data obtained by parsing the xml document
	 * @return output, a table arranged in the proper way for Filemaker to import q-values
	 * 		   based on experiment row and timepoint column
	 */

	double[][] makeSILACList(double[][] array) throws BadDataException{
		if(array == null){
			throw new BadDataException("No data with which to make the R input", path_to_wait, _error);
		}
		// output (table) dimensions = tpts x data table length -> data table has one row for each
		// experiment, and num tpts for each tpt of the experiment						
		double[][] output = new double[num_timepoints][array[0].length];
		
		
		for(int i = 0; i < array[0].length; i++){							// for each experiment (each row in Filemaker)
			for(int j = 0; j < num_timepoints; j++){						// get the p-val for each timepoint:
				double[][] sendingTable = new double[2][num_replicates];	// sendingTable = 2D array with tpt's replicate vals to be used in R
				int ind = 0;												// ind = sendingTable's current row index
				for(int k = 0; k < num_replicates; k++){					// put the replicates into the table:
					// for the silac numerator
					int lookingAtSNum = (k*num_timepoints) + (j+1) + ((num_timepoints*num_replicates)*(silac_weight_1-1));	// arithmetic explained in labelFreeCalculations
					if(!(Double.isNaN(array[lookingAtSNum][i]))){			// if there's an actual number in this part of the table
						if(log_silac){										// if the command line determined the log is to be taken
							sendingTable[0][ind] = (Math.log(array[lookingAtSNum][i])/Math.log(2));
							// see comment for this in makeLabelFreeList
						}
						else{
							sendingTable[0][ind] = array[lookingAtSNum][i];		// put the replicate's value into the table							
						}
					}
					else{													// if there's not a number
						sendingTable[0][ind] = -1;							// put placeholders into the table (TalkToR takes care of them)
					}
					// silac denominator
					int lookingAtSDen = (k*num_timepoints) + (j+1) + ((num_timepoints*num_replicates)*(silac_weight_2-1)); // repeat for silac_weight_2
					if(!(Double.isNaN(array[lookingAtSDen][i]))){
						if(log_silac){
							sendingTable[1][ind] = (Math.log(array[lookingAtSDen][i])/Math.log(2));
						}
						else{
							sendingTable[1][ind] = array[lookingAtSDen][i];							
						}
					}
					else{
						sendingTable[1][ind] = -1;
					}
					ind++;													// go to the next row of sendingTable
				}
				
			
				if(_unpaired==0){
					System.out.println("\t" + "(Getting paired t-test p-values)");
					double pVal = r_talker.pairedTtest(sendingTable, min_replicates);	// get the p-val from R using the paired t-test
					output[j][i] = pVal;
				}// put this tpt's p-val in the right place in the output array
				else if(_unpaired==1){
					double pVal = r_talker.unpairedTtest(sendingTable, min_replicates);	// get the p-val from R using the paired t-test
					output[j][i] = pVal;
				}
				
				else if(_unpaired==-1){
					double pVal = r_talker.pairedTtest(sendingTable, min_replicates);	// get the p-val from R using the paired t-test
					output[j][i] = pVal;
				}
				
				
			}
			
			
		}
		// only send actual p-values to R (no -1s), so make a list of JUST actual p-values
		java.util.ArrayList<Double> pList = new java.util.ArrayList<Double>();
		for(int j = 0; j < output[0].length; j++){
			for(int i = 0; i < output.length; i++){
				if(output[i][j] > -1){				// TalkToR returns -1s when output[i][j] for some reason shouldn't have a p-value
					pList.add(output[i][j]);		// (usu. because of not enough replicates)
				}
			}
		}
		System.out.println("\t" + "(Correcting t-test p-values for multiple tests: getting q-values)");
		String file1=file_name+"_silac_pHist.png";
		String file2=file_name+"_silac_qplot.png";
		double[] qVals = r_talker.getQValues(pList,file1,file2,BHfdr);
		// put the q-values back into output, in the spots without -1s (so the spots from which each q-value's
		//  p-value was taken, to be put into the list going to R)
		int qIndex = 0;
		for(int j = 0; j < output[0].length; j++){
			for(int i = 0; i < output.length; i++){
				if(output[i][j] > -1){
					output[i][j] = qVals[qIndex];
					qIndex++;
				}
			}
		}
		// at this point, output has a qvalue (or a placeholder-- a -1) for every timepoint (or missing timepoint)
		//	in the experiment list (so, for every possible timepoint for every Filemaker row)
		
		return output;
				
		
	}
	
	/**
	 * consolidate:
	 *    This method gets the two halves of the destined-for-Filemaker output from
	 *   labelFreeCalculations and silacCalculations and puts the two 2D arrays together. If
	 *   SILAC data was not included in the XML document input by the user, this method will 
	 *   not make the SILAC part of the total 2D array. The entire array it creates can then be 
	 *   sent on to be printed in the proper way for Filemaker.
	 *    This method can be used in general for any task involving adding two tables of identical
	 *   vertical length side by side, into one table.
	 * 	
	 * @param labelFree, the label free array constructed in labelFreeCalculations
	 * 		  SILAC, the SILAC array constructed in silacCalculations
	 * @return output, a 2D String array of the two input arrays side-by-side, consolidated into
	 * one array
	 */

	double[][] consolidate(double[][] labelFree, double[][] SILAC, double[][] data, boolean hasSilac) throws BadDataException{
		if(labelFree == null || SILAC == null){
			throw new BadDataException("No data with which to make the R input", path_to_wait, _error);
		}
		
		// output will be the 2D array that will eventually be printed for input into Filemaker
		double[][] output = new double[(labelFree.length + SILAC.length)+1][labelFree[0].length];		
		for(int i = 0; i < labelFree[0].length; i++){					// add in all the label-free values
			output[0][i] = data[0][i];									// put the counter number of the experiment in the first index				
			for(int j = 1; j < labelFree.length+1; j++){
				output[j][i] = labelFree[j-1][i];						// add into each index in this row the corresponding label-free statistic
			}
			for(int j = labelFree.length+1; j < (labelFree.length + SILAC.length + 1); j++){	// add in the SILAC values (if any)
				output[j][i] = SILAC[j-(labelFree.length+1)][i];								// (if none, this loop won't be entered)
			}
		}
		return output;
	}
	
	/**
	 * printOutputArray:
	 *   This method is for testing, to print out into a file any 2D array constructed in the code 
	 *   (as there are quite a few of them).
	 * 	
	 * @param table, any 2-dimensional array made by this code
	 * @return void
	 */

	void printOutputArray(double[][] table){			// This method is for testing, to print out any 2D array constructed in the code
		java.lang.StringBuffer buffer = new java.lang.StringBuffer();
		// Build the table
		System.out.println("table length " + table[0].length);		//144226
		for(int row = 0; row < table.length; row++){
		   for(int col = 0; col < table[0].length; col++){
			
				//	the following is to keep track of where this method is in building the array to print
				/*if(col % 500 == 0){
					System.out.println(row + ", " + col);
				}*/
				buffer.append("," + table[row][col]);					
			}
			buffer.append("\r\n");
		}
		String output = buffer.toString();
		try{
			String[] forPTF = path_to_wait.split("\\\\");
			int endLen = forPTF[forPTF.length-1].length();
			String pathToFile = path_to_wait.substring(0, path_to_wait.length()-endLen) + "consolidated.csv";
			// Create file 
			FileWriter fstream = new FileWriter(pathToFile);
			BufferedWriter out = new BufferedWriter(fstream);
			out.write(output);
			//Close the output stream
			out.close();
		}catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
	}
		
	/**
	 * printFilemakerInput:
	 *   This specifically prints output_list in order for the data to be entered into
	 *   Filemaker. It doesn't need output_list as a parameter, since output_list is an
	 *   instance variable. It does take the file name of the xml file with the original
	 *   data, in order to give the q values list that name-- this is so the q values for
	 *   a certain data set are easy to find.
	 * 	
	 * @param file, the name of the XML file from which the data to be input into Filemaker came
	 * @return void <- but it does print out a file to the same directory from which the
	 *                 xml file came
	 */
	
	void printFilemakerInput(String file, double[][] input){
		// Build the table
		java.lang.StringBuffer buffer = new java.lang.StringBuffer();
		for(int col = 0; col < input[0].length; col++){
			for(int row = 0; row < input.length; row++){
				if(row==0){					// for the counter number
					buffer.append("\t" + (int) input[row][col]);
				}
				else{						// for all the other indices (the ones with q values (or not) in them)		
					if(input[row][col] < 0){		// if this index is for a timepoint without a q value
						buffer.append("\t");		// don't print anything out at all
					}
					else{								// otherwise it's for a timepoint that has a q value
						buffer.append("\t" + input[row][col]);				// so print out that q value
					}
				}
			}
			buffer.append("\r\n");
		}		
		String output = buffer.toString();
		// Make the new filename the same as the XML filename, but a .txt file
		String name = file.replaceAll("xml", "txt").replaceAll("XML", "TXT");
		try{
			// Create file 
			FileWriter fstream = new FileWriter(name);
			BufferedWriter out = new BufferedWriter(fstream);
			out.write(output);
			//Close the output stream
			out.close();
		}catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
	}
	
	/**
	 * calcStDevs:
	 *   This method calculates the standard deviations for the measurements at each timepoint. It takes
	 *   the sums of the replicate measurements and the number of replicates summed, calculates their
	 *   average, and uses that along with the references to the other replicates in each timepoint
	 *   to calculate the standard deviation for each timepoint. If "logLF" was entered at the command line,
	 *   this method takes the logs of the label-free data (it is only called from the label-free
	 *   calculation method) to calculate the standard deviations.
	 * 	
	 * @param array, the array of data from the xml document
	 * 		  row, the row being looked at in the array of data (need this because the data is looked
	 *             at by row in makeLabelFreeList and makeSILACList, and they call this method)
	 *        sums, the list of sums of the replicate measurements at each timepoint
	 *        numSummed, the list of the number of replicates summed at each timepoint
	 *        weight, the label free weight or SILAC weight from which the replicates came for these timepoints
	 * @return stDevs, the list of standard deviations for the measurements at each timepoint
	 */

	double[] calcStDevs(double[][] array, int row, double[] sums, double[] numSummed, int weight){
		double[] avgs = new double[sums.length];
		double[] stDevs = new double[sums.length];
		for(int i = 0; i < avgs.length; i++){
			if(numSummed[i] > 0){
				avgs[i] = sums[i]/numSummed[i];
			}		// don't need an else case because make(Whichever)List takes care of the case where nothing was summed
					//	(by not making a row of the list for that result)
		}
		for(int i = 0; i < stDevs.length; i++){
			// for each timepoint
			for(int k = 0; k < num_replicates; k++){	
				int lookingAt = (k*num_timepoints) + (i+1) + ((num_timepoints*num_replicates)*(weight-1));
				if(!(Double.isNaN(array[lookingAt][row]))){			// if the value is an actual value
					double toAdd = array[lookingAt][row] - avgs[i];		// subtract the mean from the value
					if(log_lf){
						toAdd = (Math.log(array[lookingAt][row])/Math.log(2)) - avgs[i];
					}
					toAdd = (toAdd*toAdd);								// square that
					stDevs[i] += toAdd;									// get total sum
				}
			}
		}
		for(int i = 0; i < stDevs.length; i++){
			stDevs[i] = stDevs[i] / (numSummed[i] - 1);					// divide sum by n-1
			stDevs[i] = java.lang.Math.sqrt(stDevs[i]);					// take square root
		}
		return stDevs;
	}
	
	/**
	 * wait:
	 *   This method prints a "1" to a file if this program is still running, or a "0" if it has
	 *   finished and all the q values are calculated.
	 * 	
	 * @param busy, a boolean saying whether or not the program is running and still calculating
	 *              q values
	 * @return void <- but it does print out a file, "wait.txt," to whatever directory the .jar for
	 *                 this program is in
	 */

	void wait(boolean busy){
		String output;
		// Determine what to write into the file based on whether the program is busy or not
		if(busy){
			output = "1";			// if this program is running, the file will say "1"
		}
		else{
			output = "0";			// once this program is done running, the file will say "0"
		}
		// Make the filename (for the same directory that the .jar is in)
		try{
			FileWriter fstream = new FileWriter(path_to_wait);
			BufferedWriter out = new BufferedWriter(fstream);
			out.write(output);
			//Close the output stream
			out.close();
		}catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
	}
	
	void reportError(String msg){					// I guess this is for testing... @author kyu
		System.out.println("Error: " + msg);
	}
	
	/**
	 * BadDataException:
	 *   This Exception class is to catch any problems that may occur with input data, if there
	 *   wasn't enough data put in or if it's not in the right format, etc. This (obviously) should
	 *   not happen, but in case it does this Exception will give the most appropriate message for
	 *   what went wrong.
	 * 	
	 * @param msg, the message saying what the error was
	 * 		  filePath, the path to the status file
	 * 		  write, a boolean with whether or not to write the error description to a .txt file
	 * @return void <- if _error is true, it will print a .txt file with the error and stack trace
	 */

	class BadDataException extends Exception{
		public BadDataException(String msg, String filePath, boolean write){
			super(msg);
			if(write){
				this.reportException(msg, filePath);
			}
		}
		
		/**
		 * reportException:
		 *   This method prints a the Exception msg of the BadDataException and the stack trace to a file
		 *   called "error.txt," located in the same directory where the status file is
		 * 	
		 * @param msg, the Exception's message of what went wrong
		 * 		  filePath, the path to the status file
		 * @return void <- but it does print out a file, "error.txt," to whatever directory the  status
		 *                 file is in
		 */
		
		void reportException(String msg, String filePath){
			System.out.println("Error: " + msg);
			msg += "\r\n";
			msg += this.getStackTraceAsString(this);
			try{
				String[] forPTF = path_to_wait.split("\\\\");
				int endLen = forPTF[forPTF.length-1].length();
				String pathToFile = filePath.substring(0, filePath.length()-endLen) + "error.txt";
				FileWriter fstream = new FileWriter(pathToFile);
				BufferedWriter out = new BufferedWriter(fstream);
				out.write(msg);
				//Close the output stream
				out.close();
			}catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}
		}
		
		/**
		 * Gets the exception stack trace as a string.
		 * @param exception
		 * @return
		 * @author Narendra Naidu from Tech Talk blog (http://www.narendranaidu.com)
		 */
		public String getStackTraceAsString(Exception exception)
		{
			java.io.StringWriter sw = new java.io.StringWriter();
			java.io.PrintWriter pw = new java.io.PrintWriter(sw);
			pw.print(" [ ");
			pw.print(exception.getClass().getName());
			pw.print(" ] ");
			pw.print(exception.getMessage());
			exception.printStackTrace(pw);
			return sw.toString();
		}
	}
}