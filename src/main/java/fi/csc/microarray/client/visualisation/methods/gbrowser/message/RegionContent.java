package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.LinkedHashMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

/**
 * Content for given genomic region. Content is data dependent, but basically it is data parsed from tabular data. All the 
 * rows fall within the genomic region.
 *
 */
public class RegionContent implements Comparable<RegionContent> {
	
	public BpCoordRegion region;
	public LinkedHashMap<ColumnType, Object> values;

	public RegionContent(BpCoordRegion region, LinkedHashMap<ColumnType, Object> values) {
		this.region = region;
		this.values = values;
	}

	public RegionContent(BpCoordRegion region, Object concisedValue) {
		this.region = region;
		this.values = new LinkedHashMap<ColumnType, Object>();
		this.values.put(ColumnType.VALUE, concisedValue);
	}

	public int compareTo(RegionContent other) {

		int regionComparison = this.region.compareTo(other.region);

		if (regionComparison != 0) {
			return regionComparison;					
		} else {
			
			return values.toString().compareTo(other.values.toString());
		}
	}
	
	@Override
	public int hashCode() {
		return region.hashCode();
	}

	@Override
	public boolean equals(Object o) {
		RegionContent other = (RegionContent) o; 
		return region.equals(other.region) && values.toString().equals(other.values.toString());		
	}
	
	@Override
	public String toString() {
		String extra = "";
		for (Object value : values.values()) {
			extra += "\t" + value;
		}
		return region.toString(true) + extra;
	}
}
