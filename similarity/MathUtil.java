import java.lang.Float;
import java.lang.Math;
import java.util.Vector;


public class MathUtil {

	public static float nanThresh = (float)1.0;


	public MathUtil() {

	}


	public static float[][] computeSquareDist(float data[][]) {
		float[][] result = new float[data.length][data.length];

		float currSum=0;

	    for(int i=0; i < data.length; i++) {
			for(int j=i+1; j < data.length; j++) {
				currSum=0;
				for(int k=0; k < data[i].length; k++)
						currSum += (data[k][i]-data[k][j])*(data[k][i]-data[k][j]);
				result[i][j]=currSum;
			}
			System.out.println("Finished "+i);
		}

	    for(int i=0; i < data.length; i++) {
			for(int j=i+1; j < data.length; j++) {
				result[j][i]=result[i][j];
			}
		}
		return result;

	}


	public static float[][] computeInnerProducts(float data[][]) {
		float[][] result = new float[data.length][data.length];

		float currSum=0;

	    for(int i=0; i < data.length; i++) {
			for(int j=i+1; j < data.length; j++) {
				currSum=0;
				for(int k=0; k < data[i].length; k++) {
					if(! (Float.isNaN(data[i][k]) || Float.isNaN(data[j][k])) )
						currSum += data[i][k]*data[j][k];
				}
				result[i][j]=currSum;
			}
		}

	    for(int i=0; i < data.length; i++) {
			for(int j=i+1; j < data.length; j++) {
				result[j][i]=result[i][j];
			}
		}
		return result;

	}


	public static float[][] computeCosine(float data[][]) {
		float[][] result = new float[data.length][data.length];

		float[] norms = new float[data.length];
		float sum=0;
		for(int i=0; i < norms.length; i++) {
			sum=0;
			for(int j=0; j < data[i].length; j++) {
				sum += data[i][j]*data[i][j];
			}
			norms[i]=(float)Math.sqrt(sum);
		}

		float currSum=0;

	    for(int i=0; i < data.length; i++) {
			for(int j=i+1; j < data.length; j++) {
				currSum=0;
				for(int k=0; k < data[i].length; k++) {
					if(! (Float.isNaN(data[i][k]) || Float.isNaN(data[j][k])) )
						currSum += data[i][k]*data[j][k];
				}
				result[i][j]=currSum/(norms[i]*norms[j]);
			}
		}

	    for(int i=0; i < data.length; i++) {
			for(int j=i+1; j < data.length; j++) {
				result[j][i]=result[i][j];
			}
		}
		return result;

	}


	public static float[][] computeCorrelation(float data[][]) {
		float[][] result = new float[data.length][data.length];
		int[] use_flags = new int[data[0].length];
		int totNaN=0;
		float sumX=0;
		float sumY=0;
		float sumXY=0;
		float sumX2=0;
		float sumY2=0;
		int totLength;
		float den=0;


	    for(int i=0; i < data.length; i++) {
			for(int j=i+1; j < data.length; j++) {
				totNaN=0;
				sumX=0;
				sumY=0;
				sumXY=0;
				sumX2=0;
				sumY2=0;

				for(int k=0; k < data[i].length; k++) {
					if(Float.isNaN(data[i][k]) || Float.isNaN(data[j][k])) {
						use_flags[k]=0;
						totNaN++;
					}
					else {
						use_flags[k]=1;
						sumX += data[i][k];
						sumY += data[j][k];
						sumXY += data[i][k]*data[j][k];
						sumX2 += data[i][k]*data[i][k];
						sumY2 += data[j][k]*data[j][k];
					}
				}

			   totLength = data[i].length-totNaN;
		  	   if(totNaN > (nanThresh*data[i].length) || totLength < 3 ) result[i][j] = Float.NaN;
		 	   else {  // compute correlation
				   den = (float)(Math.sqrt(totLength*sumX2-sumX*sumX)*Math.sqrt(totLength*sumY2-sumY*sumY));
				   if(den < .00001) result[i][j] = Float.NaN;
				   else result[i][j] = (float)(totLength*sumXY - sumX*sumY)/den;
			   }
			}
		}

		return result;

	}

	public static float[] computeCorrelationLinearized(float data[][]) {
		int totNaN=0;
		float sumX=0;
		float sumY=0;
		float sumXY=0;
		float sumX2=0;
		float sumY2=0;
		int totLength;
		float den=0;

		int width = data.length;
		float[] result = new float[width*(width-1)/2];

		int currInd=0;

	    for(int i=0; i < data.length; i++) {
			for(int j=i+1; j < data.length; j++) {
				totNaN=0;
				sumX=0;
				sumY=0;
				sumXY=0;
				sumX2=0;
				sumY2=0;

				for(int k=0; k < data[i].length; k++) {
					if(Float.isNaN(data[i][k]) || Float.isNaN(data[j][k])) {
						totNaN++;
					}
					else {
						sumX += data[i][k];
						sumY += data[j][k];
						sumXY += data[i][k]*data[j][k];
						sumX2 += data[i][k]*data[i][k];
						sumY2 += data[j][k]*data[j][k];
					}
				}

			   totLength = data[i].length-totNaN;
		  	   if(totNaN > (nanThresh*data[i].length) || totLength < 3) result[currInd] = Float.NaN;
		 	   else {  // compute correlation
		 	   	   den = (float)(Math.sqrt(totLength*sumX2-sumX*sumX)*Math.sqrt(totLength*sumY2-sumY*sumY));
		 	   	   if(den < .00001) result[currInd] = Float.NaN;
				   else result[currInd] = (float)(totLength*sumXY - sumX*sumY)/den;
			   }
			   currInd++;
			}
		}

	  return result;
	}

/*
	public static float[][] computeZScoreCombination(Vector dataVect) {

		float[][] finalResult = new float[((float[][])dataVect.get(0)).length][((float[][])dataVect.get(0)).length];
		System.out.println("Got here");

		float[][] currResult = new float[((float[][])dataVect.get(0)).length][((float[][])dataVect.get(0)).length];

		for(int i=0; i<dataVect.size(); i++) {
			System.out.println(currResult.length + " " + currResult[0].length);
			computeCorrelation((float[][])dataVect.get(i),currResult);

			float mean=0;
			float std=0;
			for(int j=0; j < currResult.length; j++)
				for(int k=j+1; k < currResult.length; k++)
					mean += currResult[j][k];
			mean /= currResult.length*(currResult.length-1)/2;

			for(int j=0; j < currResult.length; j++)
				for(int k=j+1; k < currResult.length; k++)
					std += Math.pow((currResult[j][k]-mean),2);

			std /= currResult.length*(currResult.length-1)/2;

			for(int j=0; j < currResult.length; j++)
				for(int k=j+1; k < currResult.length; k++)
					finalResult[j][k] += (currResult[j][k]-mean)/std;
		}

		return finalResult;
	}
*/

	public static float[] computeZScoreCombination(Vector dataVect) {
		int width = ((float[][])dataVect.get(0)).length;

		float[] finalResult = new float[width*(width-1)/2];
		float[] currResult = null;


		for(int i=0; i<dataVect.size(); i++) {
			currResult = computeCorrelationLinearized((float[][])dataVect.get(i));

			float mean=0;
			float std=0;
			float tot=0;

			for(int j=0; j < currResult.length; j++)
					if(! Float.isNaN(currResult[j])) { mean += currResult[j]; tot++;}

			mean /= tot;

			for(int j=0; j < currResult.length; j++)
					if(! Float.isNaN(currResult[j])) std += Math.pow((currResult[j]-mean),2);

			std = (float)Math.sqrt(std/tot);

			for(int j=0; j < finalResult.length; j++)
					if(! Float.isNaN(currResult[j])) finalResult[j] += (currResult[j]-mean)/std;

			currResult=null;
			System.gc();
		}

		return finalResult;
	}


   public static int getIndex(int width,int i,int j) {
	   return ((2*width*i-i*i-3*i+2*j-2)/2);
   }

}