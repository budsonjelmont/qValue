import java.io.BufferedWriter;
import java.io.FileWriter;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

/**
 *	Class XMLToDataArray2:
 *		This class does the XML parsing for the program QValue, reading in an XML file and making it into 
 *		a 2D double array for class QValues.
 */

public class XMLToDataArray2 {

	private java.util.ArrayList<java.util.ArrayList<Double>> dataSList;
	private int timepoints;
	private double minPeakAr;

	public XMLToDataArray2(int tpts, double minPA){
	    timepoints = tpts;
		minPeakAr = minPA;
	}

	/**
	 * parseFilemakerData:
	 *    	This method takes an XML file (named by "file") and parses all its data to be put into a 2D
	 *    double array. It determines mathematically which kind of data it is encountering (a counter
	 *    number, a peakarea, etc) and inputs the data into the correct place in the array. In-line
	 *    comments provide a much more detailed description of the process.
	 *    	If a minimum peak area was provided by QValues, this method will put a Double.NaN into any
	 *    index in the array that would have been given a peak area lower than the minimum.
	 * 	
	 * @param file, the XML file with all the data
	 * @return dataS, the 2D data array to be used by QValues
	 */
	
	double[][] parseFilemakerData(String file){
		try {
			dataSList = new java.util.ArrayList<java.util.ArrayList<Double>>();		// this list holds lists of values, one list per experiment	
			
			SAXParserFactory factory = SAXParserFactory.newInstance();
			SAXParser saxParser = factory.newSAXParser();
			DefaultHandler handler = new DefaultHandler() {				// handler gets all specifications of what to do when seeing tags in the XML file
				boolean afterData = false;											// boolean to say if a "character" comes after "<DATA>"
				java.util.ArrayList<Double> listOfNos = new java.util.ArrayList<Double>();	// a list to hold the values from each <COL>
				int sizeOfListOfNos = 0;									// the size that the list SHOULD be if including blank <DATA> spots
																			//  (need this b/c handler won't call "characters()" for  <DATA> </DATA> spots)
				java.util.ArrayList<Double> row = new java.util.ArrayList<Double>();	// a list to hold the values of each experiment
																						//		-> this is what will go into dataSList
				public void startElement(String uri, String localName,String qName, Attributes attributes) throws SAXException {
					if (qName.equalsIgnoreCase("COL"))	{		// at <COL>, it's either a counter # or a new replicate's values
						listOfNos.clear();						// clear the list to put in the new numbers for that replicate (or the counter)
						sizeOfListOfNos = 0;					// and then reflect that the list has been cleared
					}
					if (qName.equalsIgnoreCase("DATA"))	{		// at <DATA>, there's a value after 
						afterData = true;						// mark that you need to record what "character" comes next
						sizeOfListOfNos++;						// you're SUPPOSED to add a value so increase the number
					}											// 		that's for checking purposes later: if you don't put a value,
																//		 you'll put a placeholder
				} // (don't need anything for <ROW>)

				public void endElement(String uri, String localName, String qName) throws SAXException {
					if	(qName.equalsIgnoreCase("COL"))	{				// if you find </COL>, either you just recorded counter number or
																		//  you just finished recording timepoint values for one replicate
						if(listOfNos.size()==1){						// the list will be 1 long if you recorded the counter number
							// if this is the counter number (always record)
							row.add(listOfNos.get(0));					//  add the value to the experiment list
						}
						else{											// the list will be more than 1 long if you recorded timepoint vals
							for(int i = 0; i < timepoints; i++){		// for each tpt (there may be more data placeholders in the XML file
																		//  than there are tpts -> and they'll occur after the actual tpt vals)
								if(listOfNos.get(i) >= minPeakAr){		//  add timepoint vals to the row IF they are more than the user-input
									row.add(listOfNos.get(i));			//  min peak area: if below that val, signify a measurement due to noise 
								}														
								else{									// include (or add!) timepoint placeholders (NaN) in the list		
									row.add(Double.NaN);				// -> not including XML placeholders (that's why only interating over
								}										//    listOfNos rather than the XML doc's entire <DATA> list after <COL>)
							}
						}
					}
					if	(qName.equalsIgnoreCase("DATA"))	{			// if you find </DATA>, you were supposed to add a number to the list
						if(sizeOfListOfNos > listOfNos.size()){			// if the size the list is supposed to be is bigger than it is, you 
																		//  know you should have added a value but there wasn't one in the
							listOfNos.add(Double.NaN);					//  file (looks like: <DATA></DATA>): so put a placeholder (NaN)
						}
						afterData = false;								// no longer after <DATA> (at this point it's after </DATA>)
					}													// since false, don't record any more characters until you're supposed to
					if	(qName.equalsIgnoreCase("ROW"))	{				// if you find </ROW>, you know the values from that experiment are over
						java.util.ArrayList<Double> temp = new java.util.ArrayList<Double>();	// this list is going into DataSList
						// the list "row" has everything that FileMaker spit out, but you need only as much as is in the actual experiment
						for(int i = 0; i < row.size(); i++){		// put the values in row into the list temp, but only as many as needed
							temp.add(row.get(i));
						}
						dataSList.add(temp);						// add the experiment list temp into the list of all experiments' data
						row.clear();								// get row ready for the next experiment's data						
					}
				}

				public void characters(char ch[], int start, int length) throws SAXException {
					if(afterData){									// if the character is found after <DATA>, it should be recorded
						String toRecord = new String(ch, start, length);	// this is the value to record, but it's a string
						double value = Double.parseDouble(toRecord);		// make the string into a double
						if(value==0){										// specifically for the purposes of this program, val of 0
							value = Double.NaN;								//  signifies incorrectly-recorded (or not-recorded) replicate
						}													//  (could probably take this out and set minPeakAr default to 0)
						listOfNos.add(value);								// add it to the list of values to be later added to the
					}														//  experiment's list
				}
			};
			saxParser.parse(file, handler);												// parse the file with the handler with the above specs	
			double[][] dataS = new double[dataSList.get(0).size()][dataSList.size()];	// make the 2D array of data for QValues program
			for(int j = 0; j < dataSList.size(); j++){
				for(int i = 0; i < dataSList.get(j).size(); i++){				// (each row = same length: can use dataSList.get(j).size() for all rows)
					dataS[i][j] = dataSList.get(j).get(i);						// make the data array to be returned
				}							
			}
			return dataS;

		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	void printOutputArray(double[][] table){			// This method is for testing, to print out any 2D array constructed in this code
		// Build the table
		String output = "";
		System.out.println("len " + table[0].length);
		for(int col = 0; col < table[0].length; col++){
			for(int row = 0; row < table.length; row++){
				//	the following is to keep track of where this method is in building the array to print
				//if(col % 500 == 0){
				//	System.out.println(row + ", " + col);
				//}
				output += "\t" + table[row][col];					
			}
			output += "\r\n";
		}
		try{
			String pathToFile = "C:\\Users\\qinqin\\Desktop\\test2.txt";
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


}
