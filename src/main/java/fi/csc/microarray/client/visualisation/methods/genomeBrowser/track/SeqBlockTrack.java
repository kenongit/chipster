package fi.csc.microarray.client.visualisation.methods.genomeBrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.View;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.Content;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RegionContent;

public class SeqBlockTrack extends Track{

	private Collection<RegionContent> reads = new TreeSet<RegionContent>();
	
	List<Integer> occupiedSpace = new ArrayList<Integer>();

	private Color color;

	private int RESOLUTION = 512;

	private Color[] charColors = new Color[] {
			new Color(64, 192, 64, 128), // A
			new Color(64, 64, 192, 128), // C
			new Color(128, 128, 128, 128), // G
			new Color(192, 64, 64, 128)    // T
	};

	
	private long maxBpLength;

	public SeqBlockTrack(View view, File file, Class<? extends AreaRequestHandler> handler, 
			ReadInstructions<?> readInstructions, Color color, long minBpLength, long maxBpLength)  {
		
		super(view, file, handler, readInstructions);		
		this.color = color;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
		
		
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		
//		Collection<Region> toBeRemoved = new ArrayList<Region>();
		
		occupiedSpace.clear();

		if(reads != null){

			
			
			Iterator<RegionContent> iter = reads.iterator();
			while(iter.hasNext()){

				RegionContent read = iter.next();
				
				Object valueObj = read.values.get(Content.VALUE);
				
				
				if(!read.region.intercepts(getView().getBpRegion())){
					
					iter.remove();
					continue;
				}
//				
//				if(valueObj == null){
//					
//					System.out.println("null value obj");
//					drawables.add(createDrawable(read.region.start, read.region.end, Color.red));
//				} else
				{
					
//					int limited = Math.min((int)(Math.log((Float) valueObj * 10)*50), 255);					
//					limited = Math.max(0, limited);
//					
//					//System.out.println(limited);
//					
//					Color c = new Color(color.getRed(), color.getGreen(), color.getBlue(), limited);
					Color c = Color.gray;
					int height = 10; //(int) (getView().getTrackHeight() / 2 * limited / 255f);
					//System.out.println(height);
					
					long startBp = read.region.start;
					long endBp = read.region.end;			
					
					
					String seq = ((String)read.values.get(Content.SEQUENCE));
					
					
					if(seq != null){
						seq = seq.trim();						
					} else {
						int seqLength = (int) (endBp - startBp + 1);
						if(seqLength > 128){
							
							//TODO what happened?
							seqLength = 0;
						}
						seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
							.substring(0, seqLength);
					}
					
//					if(((String)read.values.get(Content.STRAND)).trim().equals("R")){
//						startBp -= seq.length();
//						endBp -= seq.length();
//					}
					
					Rectangle rect = new Rectangle();

					rect.x = getView().bpToTrack(startBp);
					rect.width = getView().bpToTrack(endBp) - rect.x;
					
					//System.out.println(startBp + ", " + endBp + ", " + seq + ", " + rect.width);
					

					int i = 0;

					while(occupiedSpace.size() > i && occupiedSpace.get(i) > rect.x + 1){
						i++;
					}

					int end = rect.x + rect.width;

					if(occupiedSpace.size() > i){
						occupiedSpace.set(i, end);
					} else {
						occupiedSpace.add(end);
					}

					rect.y = (int)(getView().getTrackHeight() - ((i + 1) * (height + 2)));
					rect.height = height;	

					//String seq = ((String)read.values.get(Content.SEQUENCE)).trim();
					
					String strand = ((String)read.values.get(Content.STRAND));
					
					if(strand != null){
						strand = strand.trim();						
					} else {
						strand = "F";
					}
					
					if(strand.equals("R")){
						StringBuffer buf = new StringBuffer(seq).reverse();
						
//						for(int j = 0; j < seq.length(); j++){
//							switch(buf.charAt(j)){
//							case 'A': buf.setCharAt(j, 'T'); break;
//							case 'C': buf.setCharAt(j, 'G'); break;
//							case 'G': buf.setCharAt(j, 'C'); break;							
//							case 'T': buf.setCharAt(j, 'A'); break;
//							}
//						}											
						
						seq = buf.toString();
						
						//seq = (String)read.values.get(Content.QUALITY);
						
					}
					
					if(rect.width < seq.length()){
						
						drawables.add(new RectDrawable(rect, c, null));
						
					} else { //Enough space to show the actual sequence
						
						final int CHAR_WIDTH = 7;
						
						float x = rect.x;							
						float increment = (rect.width)/((float)seq.length());
						
												
						for(int j = 0; j < seq.length(); j++){
							
							char letter = seq.charAt(j);
							
							if(rect.width > seq.length() * CHAR_WIDTH){
								
								drawables.add(new TextDrawable((int)x, rect.y + 8, 
										"" + letter, Color.black));
								
//								drawables.add(new TextDrawable((int)x, rect.y - 8, 
//										"" + letter, Color.black));
							
							}
							
							Color bg = Color.white;
							if(letter == 'A'){
								bg = charColors[0];
							} else if (letter == 'C'){
								bg = charColors[1];
							} else if (letter == 'G'){
								bg = charColors[2];	
							} else if (letter == 'T'){
								bg = charColors[3];
							}
							
							drawables.add(new RectDrawable((int)x, rect.y - 2, (int)increment , 10,bg,null));
							
							x += increment;
						}
					}			
				}				
			}
		}
				
		return drawables;
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {		

//		if (areaResult.content instanceof List) {
//			List<List<Object>> reads = (List<List<Object>>) areaResult.content;
//
//			for(List<Object> obj: reads){
//
//				//TODO make separate track for database genes and read and remove this hack 
//				long start = (Long)obj.get(areaResult.fileDef.indexOf(Content.BP_START)); 
//				long end;
//				if(areaResult.fileDef.indexOf(Content.BP_END) >= 0){
//					end = (Long)obj.get(areaResult.fileDef.indexOf(Content.BP_END)); 
//				} else {
//					end = start + ((String)obj.get(areaResult.fileDef.indexOf(Content.SEQUENCE))).length();
//				}
//
//				Region reg = new Region(start, end);
//
//				this.reads.add(new RegionValue<Float>(reg, 100f));
//			}			
//			
//		} else 

		if(areaResult.status.concise == isConcised() && 
				isForForwardStrand(areaResult) != isReversed()){

			this.reads.add(areaResult.content);			

			getView().redraw();
		}

		//this.reads.addAll(result.collection);
	}

	private boolean wasLastConsied = true;

	private long minBpLength;
	
	@Override
	public void updateData(){

		if(wasLastConsied != isConcised()){
			reads.clear();
			wasLastConsied = isConcised();
		}
		super.updateData();
	}
	
	@Override
	public int getMaxHeight(){
		if(getView().getBpRegion().getLength() > minBpLength && 
				getView().getBpRegion().getLength() <= maxBpLength){
			
			return super.getMaxHeight();
		} else {
			return 0;
		}
	}

	@Override
	public Collection<Content> getDefaultContents() {
		return Arrays.asList(new Content[] {
				Content.SEQUENCE, Content.STRAND, Content.QUALITY
				}); 
	}

	@Override
	public boolean isConcised() {
		return getView().getBpRegion().getLength() > 1*1024*1024;
		//return reads.size() > RESOLUTION;
	}
}