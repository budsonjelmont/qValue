import java.io.BufferedWriter;
import java.io.FileWriter;

import org.rosuda.JRI.*;

/**
 *	Class TalkToR:
 *		This class does the Java-to-R interfacing for the program QValue.
 * 
 */

public class TalkToR {

	private Rengine _re;
	
	

	/**
	 * getQValues:
	 *    This method takes a list made in QValues, of p-values, (calculated by any statistical test)
	 *   and manually makes a list in R of this data. It tells R to load the library of the "qvalue"
	 *   program, and use the qvalue function of this program to calculate the q-values from the 
	 *   p-values. It then converts the R list of q values into a Java list, (double array) to return.
	 *   
	 *   Note: qvalue must already be installed on the machine in order for this method to work
	 * 	
	 * @param pVals, the list made in QValues for R to calculate the q-values
	 * @return output, the list of q values in the correct order to be input into the 2D array
	 * destined for Filemaker
	 */
	double[] getQValues(java.util.ArrayList<Double> pVals, String phist, String qplot, int BHfdr){
		if(_re!=null){												// if there's an REngine to work with (will only be null if R failed to load)
			try{
				String toEval = "pvals <- array(0:0, dim=" + pVals.size() + ")";					// makes a list in R
				_re.eval(toEval,false);			// now you have a list full of zeros with enough places for each experiment's p-value
				java.util.Stack<String> commandStack = new java.util.Stack<String>();	// stack so R commands can run without overflowing
				for(int i = 0; i < pVals.size(); i++){
					toEval = "pvals[" + (i+1) + "] <- " + pVals.get(i);					// put each p-value in the list into the list made in R
					commandStack.push(toEval);											// fill the stack
				}
				while(!commandStack.isEmpty()){							// use the created stack like the memory stack and make the table in R
					String command = commandStack.pop();
					_re.eval(command);
				}
				// "pvals" is the p value list : get the q values using the R function qvalue
				_re.eval("library(qvalue)");							// load the qvalue library and use it to calculate q values
				
				//If BHfdr was specified by the user, set lambda to 0 so that pi0 becomes 1 (this is the same as controlling FDR with the Benjamini-Hochberg method)
				if(BHfdr==1){
					System.out.println("Controlling FDR with the Benjamini-Hochberg method");
					_re.eval("qobj <- qvalue(pvals, lambda=0)");
				} else {
					_re.eval("qobj <- qvalue(pvals)");
				}
				
				REXP qvals = _re.eval("qvals <- qobj$qvalue");		// qvals is the list in R of all the q values
				if(qvals==null){									// if qvals is null, it's because pi0 couldn't be estimated with the default lambda parameters
					_re.eval("qobj <- qvalue(pvals, pi0.method=\"bootstrap\")");	// first try to use the bootstrap method of estimating pi0
					 qvals = _re.eval("qvals <- qobj$qvalue");
					
				}
				// if that doesn't work, the while loop below will be entered: try changing the upper bound of lambda list for estimating pi0
				double upper = .95;
				while(qvals==null){															// loop through until you get a real pi0
					_re.eval("qobj <- qvalue(pvals, lambda=seq(0, " + upper + ", .01))");
					qvals = _re.eval("qvals <- qobj$qvalue");								// decrease the upper bound each time to try again if
					
					upper -= .01;															//  pi0 isn't estimated and thus a list isn't returned
					if(upper < .01){														// if upper has been changed to 0, don't check again
						break;																//  because the lambda list will just be from 0 to 0
					}																		// break out of the while loop and continue
				}
				// if qvals is still null, while loop below will be entered: change the lower bound of lambda list for estimating pi0 until you get a real pi0
				double lower = 0;
				while(qvals==null){															
					_re.eval("qobj <- qvalue(pvals, lambda=seq(" + lower + ", .99, .01))"); // increase the upper bound each time to try again if
					qvals = _re.eval("qvals <- qobj$qvalue");								//  pi0 isn't estimated and thus a list isn't returned
					
					lower += .01;															// if upper has been changed to .99, don't check again
					if(lower > .98){														//  because lambda list will just be from .99 to .99
						break;																// if that happens, break out of while loop
					}
				}
				// in the end if qvals is still null - if NONE of that worked - just use pi0 = 1 (same as the stepwise Benjamini and Hochberg method)
				if(qvals==null){
					_re.eval("qobj <- qvalue(pvals, lambda=0)");							// when lambda = 0, pi0 automatically set = 1
					qvals = _re.eval("qvals <- qobj$qvalue");
					
					
				}
				_re.assign("phist",phist);
				_re.eval("png(filename=phist)");
				_re.eval("hist(pvals,breaks=100)");
				_re.eval("dev.off()");
				_re.assign("qplot",qplot);
				_re.eval("png(filename=qplot)");
				_re.eval("qplot(qobj)");
				_re.eval("dev.off()");
				
				
				_re.eval("data<-cbind(pvals,qvals)");
				
				//Write csv of p&q values for troubleshooting
				//String name="C:/Users/superuser.DCMASTER.001/Desktop/p_q.csv";
				//_re.assign("sourceName", name);
				//_re.eval("write.csv(data, file=sourceName)");
				
				double[] output = qvals.asDoubleArray();					// make the "qvals" list in R into an array of doubles			
				
				return output;
			} catch (Exception e) {						// if the REngine exists but there's some exception in its use
				System.out.println("EX:"+e);
				e.printStackTrace();
			}
		}
		else{											// if the REngine is null (so, R couldn't get started)
			System.out.println("Cannot load R");			
		}
		return null;									// so if R doesn't open, it won't return any array
	}
	
	/**
	 *  This method takes the table made in QValues, of the sums of data, number of replicates summed,
	 *   and their standard deviations, and manually makes a table in R of this data. It tells R to
	 *   calculate p values using this table, using the small-sample test for the Student's t distribution.
	 *   It then converts the R list of p values into a Java list, (double array) to return.
	 *   
	 *   @param table, the table made in QValues for R to perform the small-sample t-test and calculate
	 * 	the p-value from this test
	 */

	double[] getSmallSamplePVals(double[][] tableForR){
		if(_re!=null){
			try{
				int numRows = tableForR[0].length;
				// make the table in R
				String toEval = "table <- array(0:0, dim=c(" + numRows + ", 6))";
				_re.eval(toEval,false);			// now you have a table full of zeros
				// put stuff in the table
				// make a stack so the R commands can run without overflowing
				java.util.Stack<String> commandStack = new java.util.Stack<String>();
				for(int i = 0; i < numRows; i++){		// for each row of the input table
					for(int j = 0; j < 6; j++){			// for each column -- there's no way it will be > or < than 6
						toEval = "table[" + (i+1) + "," + (j+1) + "] <- ";			// one-based indexing
						double num = tableForR[j][i];
						toEval += num;
						commandStack.push(toEval);
						//System.out.println(toEval);
					}
				}
				// use the created stack like the memory stack and make the table in R
				while(!commandStack.isEmpty()){
					String command = commandStack.pop();
					_re.eval(command);
				}
				// now make the variables you need to make the p value list
				_re.eval("X <- table[,1]",false);
				_re.eval("Y <- table[,2]",false);
				_re.eval("n <- table[,3]",false);
				_re.eval("m <- table[,4]",false);
				_re.eval("s.X <- table[,5]",false);
				_re.eval("s.Y <- table[,6]",false);
				// do the calculations to make the p value list
				// use small-sample test (the T statistic and the Student's t distribution)
				_re.eval("Sp <- sqrt((((n-1)*(s.X^2))+((m-1)*(s.Y^2)))/(n+m-2))");
				_re.eval("T <- ((X/n)-(Y/m))/(Sp*sqrt((1/n)+(1/m)))");
				_re.eval("df = n+m-2");
				REXP pvals = _re.eval("pvals <- 2*(pt(abs(T),df,lower.tail=FALSE))");
				
				double[] output = pvals.asDoubleArray();
				// now output is the array of all the p values
				_re.end();
				return output;
			} catch (Exception e) {
				System.out.println("EX:"+e);
				e.printStackTrace();
			}
		}
		else{
			System.out.println("Cannot load R");			
		}
		_re.end();
		return null;
	}
	
	/**
	 * pairedTtest:
	 * 	  This method takes the table made in QValues, of the two different SILAC weight replicates
	 *   for a timepoint-- column one with the replicates from one SILAC weight (the user-input 
	 *   silac_numerator) for a timepoint, and column two with the replicates from the other SILAC
	 *   weight (the user-input silac_denominator) for that timepoint-- and uses the table to 
	 *   count the number of replicate pairs for the timepoint; if this number is appropriate for
	 *   obtaining meaningful statistics, (see in-line comments) this method calls makeTable() to 
	 *   send commands to R to make the table. Then pairedTtest() interfaces with R to perform R's
	 *   t-test with the "paired" option =TRUE, getting R's return object for this paired t-test 
	 *   and interpreting it for Java, in order to extract the p-value.
	 *    If the p-value cannot be calculated for some reason (such reasons for this can be found 
	 *   in the in-line comments of this method) a -1 will be returned to QValues as a placeholder
	 *   for the p-value instead of an actual p-value.
	 * 	
	 * @param table, the table made in QValues for R to perform the paired t-test and calculate
	 * the p-value from this test
	 * 		  minReplicates, the minimum number of replicates (input by the user in the command 
	 * line or at default value of 2) needed in order to accept as meaningful the statistical output
	 * of the test -> in the case of this method, minReplicates is the minimum number of replicate
	 * PAIRS needed
	 * @return p, the p-value denoting the likelihood for a change between two different SILAC
	 * weights for a given timepoint
	 */

	double pairedTtest(double[][] table, int minReplicates){
		if(_re!=null){									// if there's an REngine to work with (will only be null if R failed to load)
			try{
				// R's paired t test will work as long as the timepoint has at least two pairs of SILAC replicates
				// As long as there are two SILAC_1 replicates that are PAIRED with SILAC_2 replicates, t.test(x,y,paired=TRUE) will run
				// But: go with whatever was set as minimum replicates needed per timepoint, as the minimum number of pairs needed (default = 2)
				
				int numPairs = 0;
				for(int i = 0; i < table[0].length; i++){
					if(table[0][i] > -1 && table[1][i] > -1){
						numPairs++;
					}
				}
				if(numPairs >= minReplicates){	// if there are enough replicates for the paired t test
												// (if # replicates >= minimum input # of replicates w/ both silac_1 and silac_2 measurements)
					this.makeTable(table, numPairs, 2);			// table for t.test has 2 cols: 
																//  1 for silac_1 replicates for this tpt, 1 for silac_2 reps for it
					// t.test(silac_1_reps, silac_2_reps, paired=TRUE) gives the return object of the paired t test,
					//  which has the resulting p-value in it (for the particular timepoint being tested)
					
					REXP resT = _re.eval("resT <- t.test(table[,1], table[,2], paired=TRUE)");		// in QValues, the table to be sent
																									//  to R is created with tpts first
																									//  and vals second, so in R table[,2]
																									//  is vals and table[,1] is tpt #				
					if(resT != null){

						// resT.asList() gives the return object's list-of-lists
						//  -> = a list-of-lists = <t value "list">,<degrees of freedom "list">,<p val "list">,<confidence interval list>,<mean of differences list>
						//  ( -> for completeness of description: list-of-lists includes also more "lists" of input parameters, after those above)
						// resT.asList().at(2) gives the list with the p-val
						//  -> list with p-val = <p val "list">
						// can be accessed as double array, (easier for Java than dealing with JRI's RList class):  resT.asList().at(2).asDoubleArray()
						//  -> list now as double array with p-val = [p val "list"]
						// resT.asList().at(2).asDoubleArray()[0] gives the first (and only) index of the p-val list, which is the p-val
						double p = resT.asList().at(2).asDoubleArray()[0];//IMPORTANT: .at(1) IS DEGREES OF FREEDOM, .at(2) IS THE P_VAL
						if(!(Double.isNaN(p))){	
							return p;
						}
						else{			// p will be NaN if for ALL replicate pairs for that experiment, silac-1 rep == silac-2 rep
							return -1;	//  (as in, if: silac-1 rep 1 = silac-2 rep 1, silac-1 rep 2 = silac-2 rep 2, etc: ALL rep pairs equal)
						}				//   return a placeholder, since no actual p-val returned by R's paired t-test	
					}
					else{			// resT will be null if the "data are essentially constant" error comes up for this table in R
						return -1;	// this means that the p-value is meaningless, so return a placeholder
					}
				}
				else{			// if there wasn't enough data to get a meaningful answer
					return -1;	//  (= if number of pairs of reps < minimum reps accepted for statistical calculations)
				}
			} catch (Exception e) {						// if the REngine exists but there's some exception in its use
				System.out.println("EX:"+e);
				e.printStackTrace();
			}
		}
		else{											// if the REngine is null (so, R couldn't get started)
			System.out.println("Error: Cannot load R");		
		}
		return -1;			// this will only be reached if R didn't load, so R not having loaded acts like no meaningful p-val returned
	}
	
	
	
	/**
	 * unpairedTtest:
	 * 	  This method takes the table made in QValues, of the two different SILAC weight replicates
	 *   for a timepoint-- column one with the replicates from one SILAC weight (the user-input 
	 *   silac_numerator) for a timepoint, and column two with the replicates from the other SILAC
	 *   weight (the user-input silac_denominator) for that timepoint-- and uses the table to 
	 *   count the number of replicate pairs for the timepoint; if this number is appropriate for
	 *   obtaining meaningful statistics, (see in-line comments) this method calls makeTable() to 
	 *   send commands to R to make the table. Then pairedTtest() interfaces with R to perform R's
	 *   t-test with the "paired" option =TRUE, getting R's return object for this paired t-test 
	 *   and interpreting it for Java, in order to extract the p-value.
	 *    If the p-value cannot be calculated for some reason (such reasons for this can be found 
	 *   in the in-line comments of this method) a -1 will be returned to QValues as a placeholder
	 *   for the p-value instead of an actual p-value.
	 * 	
	 * @param table, the table made in QValues for R to perform the paired t-test and calculate
	 * the p-value from this test
	 * 		  minReplicates, the minimum number of replicates (input by the user in the command 
	 * line or at default value of 2) needed in order to accept as meaningful the statistical output
	 * of the test -> in the case of this method, minReplicates is the minimum number of replicate
	 * PAIRS needed
	 * @return p, the p-value denoting the likelihood for a change between two different SILAC
	 * weights for a given timepoint
	 */

	double unpairedTtest(double[][] table, int minReplicates){
		if(_re!=null){									// if there's an REngine to work with (will only be null if R failed to load)
			try{
				// R's paired t test will work as long as the timepoint has at least two pairs of SILAC replicates
				// As long as there are two SILAC_1 replicates that are PAIRED with SILAC_2 replicates, t.test(x,y,paired=TRUE) will run
				// But: go with whatever was set as minimum replicates needed per timepoint, as the minimum number of pairs needed (default = 2)
				
				
				int numSummed1=0;
				int numSummed2=0;
				double sums1=0;
				double sums2=0;
				double stDevs1=0;
				double stDevs2=0;
				double mean1=0;
				double mean2=0;
				
				for (int i=0; i<table[0].length;i++){
					if(table[0][i]>-1){
						numSummed1++;
						sums1=sums1+table[0][i];						
					}				
				}
				
				for (int i=0; i<table[1].length;i++){
					if(table[1][i]>-1){
						numSummed2++;
						sums2=sums2+table[1][i];						
					}				
				}
				
				if(sums1>=minReplicates&&sums2>minReplicates ){
				mean1=sums1/numSummed1;
				mean2=sums2/numSummed2;
				
				for (int i=0; i<table[0].length;i++){
					if(table[0][i]>-1){
						double toAdd=table[0][i]-mean1;
						toAdd=(toAdd*toAdd);
						stDevs1 += toAdd;						
					}				
				}
				
				for (int i=0; i<table[0].length;i++){
					if(table[1][i]>-1){
						double toAdd=table[1][i]-mean2;
						toAdd=(toAdd*toAdd);
						stDevs2 += toAdd;						
					}				
				}
				
				stDevs1=stDevs1/(numSummed1-1);
				stDevs1 = java.lang.Math.sqrt(stDevs1);
				
				stDevs2=stDevs2/(numSummed2-1);
				stDevs2 =java.lang.Math.sqrt(stDevs2);
				
				_re.eval("X <-"+sums1+"");
				_re.eval("Y <- "+sums2+"");				
				_re.eval("n <- "+numSummed1+"");
				_re.eval("m <- "+numSummed2+"");
				_re.eval("s.X <- "+stDevs1+"");
				_re.eval("s.Y <- "+stDevs2+"");
				
				_re.eval("Sp <- sqrt((((n-1)*(s.X^2))+((m-1)*(s.Y^2)))/(n+m-2))");
				_re.eval("T <- ((X/n)-(Y/m))/(Sp*sqrt((1/n)+(1/m)))");
				_re.eval("df = n+m-2");
				REXP pvals = _re.eval("pvals <- 2*(pt(abs(T),df,lower.tail=FALSE))");	
				//REXP pvals=_re.eval("pvals <- 2*(1 - pt(abs(T),df))",false);
				if(pvals!=null){
					double p=pvals.asDouble();
					
					if(!(Double.isNaN(p))){
						_re.end();
					    return p;
					}
					else{
						return -1;
					}
				}
				else{
					return -1;
				}

			
				}
				else {
					return -1;  //not enough replicates
				}
			
			}catch (Exception e) {
				System.out.println("EX:"+e);
				e.printStackTrace();
			}
			
			
		}		
		else{
	     	 System.out.println("can not load R");
		}
		
		return  -1;
	}
	
	
	
	
	/**
	 * makeTable:
	 * 
	 * 	This method takes a table that was made in the form of a 2D double array (instantiated in 
	 * class QValues) and tells the instance of R that is running (with which this program is interfacing)
	 * to make this table in R. It is held in the variable "table" in R.
	 * 	Each time a table needs to be made, it will be made in this method in R. Each table made will
	 * overwrite the previous table, and EACH time instantiated, the size of "table" is specified (making
	 * a completely new instance) so there will never be any overlap between "table"s. Thus every time a
	 * table is made, all appropriate calculations must be performed on it immediately.
	 * 	This method checks each index and makes sure the index is not -1. If an index holds -1, makeTable()
	 * assumes this is a placeholder from QValues and it will ignore that entire row of the table. 
	 * 
	 * @param table, the table sent from QValues to be made into a table in R
	 * @param rows, the number of rows in the table sent from QValues
	 * @param cols, the number of columns to be used: usually depends on statistical test to be calculated,
	 * 				this param determined in the method calling makeTable()
	 * @return void --> makes a table in the instance of R that is running
	 */
	
	void makeTable(double[][] table, int rows, int cols){
		double[][] tableForR = new double[cols][rows];	// as many cols as were input
		
		
		
		// fill the table with all the actual data
		int row = 0;
		// make the table for everything without a -1 in either column --> (a -1 in either or both columns represents placeholders
		for(int i = 0; i < table[0].length; i++){						//  OR a pair that isn't truly a pair because one val is missing)
			if(table[0][i] > -1 && table[1][i] > -1){
				tableForR[0][row] = table[0][i];						// put the values of the input table into the table to be sent to R
				tableForR[1][row] = table[1][i];
				row++;
			}// now the table is ready to be put into R
		}
		String toEval = "table <- array(0:0, dim=c(" + rows + ", " + cols + "))";	// table constructed w/ correct number of rows (one for
		_re.eval(toEval,false);														//  each replicate) and correct number of cols (however
																					//  many necessary for desired statistical test) 
																					//  -> makes a table full of zeros in R
		// put the actual data into the table in R
		java.util.Stack<String> commandStack = new java.util.Stack<String>();	// stack so R commands run w/out stack overflow
		// fill the stack
		for(int i = 0; i < tableForR[0].length; i++){			// for each row of the input table
			for(int j = 0; j < cols; j++){							// for each column
				toEval = "table[" + (i+1) + "," + (j+1) + "] <- ";		// i+1, j+1 b/c R uses one-based indexing
				double num = tableForR[j][i];
				toEval += num;
				commandStack.push(toEval);
			}
		}
		while(!commandStack.isEmpty()){					// use the created stack like the memory stack and make the table in R
			String command = commandStack.pop();			// using it like this avoids stack overflow
			_re.eval(command);
		}												// now the table is constructed in R and called "table"
	}

	/**
	 * openR:
	 *   This method checks that the R version installed is the correct version, and that R
	 *   can be loaded and the Rengine used for this program can be made. If not, it prints out
	 *   an error message.
	 * 	
	 * @param none
	 * @return re, an Rengine to be used to communicate with R and do calculations therein
	 */

	Rengine openR(){
		if (!Rengine.versionCheck()) {
			System.err.println("** Version mismatch - Java files don't match library version.");
			System.exit(1);
		}
		// Creating the R Engine	(this is code from the "rtest" file of the JRI package)
		String[] r_args = null;
		//		Rengine re=new Rengine(r_args, false, new TextConsole());
		_re=new Rengine(r_args, false, new TextConsole());
		// the engine creates R is a new thread, so we should wait until it's ready
		if (!_re.waitForR()) {
			System.out.println("Cannot load R");
			return null;
		}
		return _re;
	}

	/**
	 * writeToFile:
	 *   This method is for testing, to print out into a file any String constructed in the code.
	 *   (Typically this String would be one constructed to show a list or a 2D array.)
	 * 	
	 * @param output, any String made by this code
	 * @return void
	 */

	void writeToFile(String output){
		try{
			// Create file 
			FileWriter fstream = new FileWriter("test.txt");
			BufferedWriter out = new BufferedWriter(fstream);
			out.write(output);
			//Close the output stream
			out.close();
		}catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}
	}

}

/**
 * class TextConsole:
 *   The class TextConsole is required to create an R Engine that is ready to be run; for the
 *   purposes of this code, TextConsole doesn't need any methods, as the R Engine used here never
 *   needs to ensure that anything is printed out for a user to see. 	
 */

class TextConsole implements RMainLoopCallbacks
{
	public void rWriteConsole(Rengine re, String text, int oType) {
	}

	public void rBusy(Rengine re, int which) {
	}

	public String rReadConsole(Rengine re, String prompt, int addToHistory) {
		return null;
	}

	public void rShowMessage(Rengine re, String message) {
	}

	public String rChooseFile(Rengine re, int newFile) {
		String res=null;
		return res;
	}

	public void   rFlushConsole (Rengine re) {
	}

	public void   rLoadHistory  (Rengine re, String filename) {
	}			

	public void   rSaveHistory  (Rengine re, String filename) {
	}			
}