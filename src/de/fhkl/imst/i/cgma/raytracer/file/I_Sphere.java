package de.fhkl.imst.i.cgma.raytracer.file;
 
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class I_Sphere extends RT_Object {
	public float[] material;
	public int materialN;
	public float[] center;
	public float radius;

	@Override
	public String getHeader() {
		return "IMPLICIT_SPHERE";
	}

	@Override
	protected void readContent(LineNumberReader br) throws IOException {
		Pattern pMaterial = Pattern.compile(materialRegex);
		Pattern pParameter = Pattern.compile(paramterRegex);
		material = new float[9];
		center = new float[3];
		
		// Material lesen
		Matcher matcher = pMaterial.matcher(readLine(br).trim());
		if(!matcher.matches())
			throw new IOException("Ungültiges Dateiformat! " + br.getLineNumber());
		for(int i = 0; i < 9; ++i)
			material[i] = Float.parseFloat(matcher.group(i+1));
		materialN = Integer.parseInt(matcher.group(10));
		// Parameter lesen
		matcher = pParameter.matcher(readLine(br).trim());
		if(!matcher.matches())
			throw new IOException("Ungültiges Dateiformat! " + br.getLineNumber());
		for(int i = 0; i < 3; ++i)
			center[i] = Float.parseFloat(matcher.group(i+1));
		radius = Float.parseFloat(matcher.group(4));
		calcBoundingBox();
	}
	
	@Override
	public void calcBoundingBox() {
                
                this.min[0] = this.center[0] - this.radius;
                this.min[1] = this.center[1] - this.radius;
                this.min[2] = this.center[2] - this.radius;
                
                this.max[0] = this.center[0] + this.radius;
                this.max[1] = this.center[1] + this.radius;
                this.max[2] = this.center[2] + this.radius;
	}

	private static final String materialRegex =
			"(\\-?[0-9]+\\.[0-9]+) +(\\-?[0-9]+\\.[0-9]+) +(\\-?[0-9]+\\.[0-9]+) +(\\-?[0-9]+\\.[0-9]+) +(\\-?[0-9]+\\.[0-9]+) +(\\-?[0-9]+\\.[0-9]+) +(\\-?[0-9]+\\.[0-9]+) +(\\-?[0-9]+\\.[0-9]+) +(\\-?[0-9]+\\.[0-9]+) +([0-9]+)";
	private static final String paramterRegex =
			"(\\-?[0-9]+\\.[0-9]+) +(\\-?[0-9]+\\.[0-9]+) +(\\-?[0-9]+\\.[0-9]+) +(\\-?[0-9]+\\.[0-9]+)";
}
