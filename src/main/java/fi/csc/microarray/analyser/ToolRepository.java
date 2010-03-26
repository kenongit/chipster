package fi.csc.microarray.analyser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.messaging.message.DescriptionMessage;
import fi.csc.microarray.messaging.message.DescriptionMessage.Category;
import fi.csc.microarray.util.XmlUtil;


/**
 * The ToolRepository manages analysis tools and runtimes. Analysis tools
 * are visible to AnalyserServer as AnalysisDescriptions, where as runtimes
 * are not visible outside of this class.
 * 
 * After initialization, any access to descriptions maps should be use 
 * synchronized using this, since reading descriptions may cause a
 * description to be updated.
 *  
 * @author hupponen
 *
 */
public class ToolRepository {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(ToolRepository.class);
	
	private LinkedHashMap<String, AnalysisDescription> descriptions = new LinkedHashMap<String, AnalysisDescription>(); 
	private LinkedHashMap<String, AnalysisDescription> supportedDescriptions = new LinkedHashMap<String, AnalysisDescription>();
	private LinkedHashMap<String, AnalysisDescription> visibleDescriptions = new LinkedHashMap<String, AnalysisDescription>();
	
	private HashMap<String, ToolRuntime> runtimes = new HashMap<String, ToolRuntime>();
	
	private DescriptionMessage descriptionMessage = null;
	
	
	/**
	 * 
	 * @param the root workDir for the jobs of the computing service
	 * @throws Exception
	 */
	public ToolRepository(File workDir) throws Exception {
		loadRuntimes(workDir);
		loadModules();
	}
	
	
	
	
	public synchronized AnalysisDescription getDescription(String fullName) throws AnalysisException {
		AnalysisDescription desc; 

		// get the description
		desc = descriptions.get(fullName);

		// check if description needs to be updated
		if (desc != null && !desc.isUptodate()) {
			updateDescription(desc);
		}
		
		// return the possibly updated description
		return descriptions.get(fullName); 
	}


	/** 
	 * @return DescriptionMessage that can be sent to the client.
	 */
	public DescriptionMessage getDescriptionMessage() {
	    return descriptionMessage;
	}


	
	/**
	 * 
	 * TODO: remove this and use getDescriptionMessage instead
	 * 
	 * Returns one huge SADL block that contains all loaded analysis 
	 * descriptions.
	 * @return huge block
	 * @throws AnalysisException 
	 */
	public synchronized StringBuffer serialiseAsStringBuffer() throws AnalysisException {
		StringBuffer buf = new StringBuffer();

		// find descs that need to be updated (custom script available)
		List<AnalysisDescription> descsToBeUpdated = new LinkedList<AnalysisDescription>();
		for (AnalysisDescription description : visibleDescriptions.values()) {
			if (!description.isUptodate()) {
				descsToBeUpdated.add(description);
			}
		}
		
		// update (can't update in the previous loop, would cause concurrent modification)
		for (AnalysisDescription description: descsToBeUpdated) {
			updateDescription(description);
		}
		
		// get the descriptions
		for (AnalysisDescription description: visibleDescriptions.values()) {
			buf.append(description.getSADL() + " ");
		}
		return buf;
	}

	public synchronized boolean supports(String fullName) {
		return supportedDescriptions.containsKey(fullName);
	}

	
	private void updateDescription(AnalysisDescription desc) throws AnalysisException {
		AnalysisDescription newDescription = desc.getHandler().handle(desc.getSourceResourceName());
		if (newDescription != null) {
			newDescription.setUpdatedSinceStartup();
			
			// name (id) of the tool has not changed
			if (desc.getFullName().equals(newDescription.getFullName())) {
				
				// replace the old description with the same name
				descriptions.put(newDescription.getFullName(), newDescription);
				if (supportedDescriptions.containsKey(desc.getFullName())) {
					supportedDescriptions.put(newDescription.getFullName(), newDescription);
				}
				if (visibleDescriptions.containsKey(desc.getFullName())) {
					visibleDescriptions.put(newDescription.getFullName(), newDescription);
				}
			} 

			// name (id) of the tool has changed
			else {
				logger.warn("name of the tool has changed after loading from custom-scripts, keeping both old and new");
				if (descriptions.containsKey(newDescription.getFullName())){
					logger.warn("descriptions already contains a tool with the new name, ignoring custom-scripts");
					return;
				} 
				// add the tool with the new name
				descriptions.put(newDescription.getFullName(), newDescription);
				if (supportedDescriptions.containsKey(desc.getFullName())) {
					supportedDescriptions.put(newDescription.getFullName(), newDescription);
				}
				if (visibleDescriptions.containsKey(desc.getFullName())) {
					visibleDescriptions.put(newDescription.getFullName(), newDescription);
				}
			}
		}
	}

	
	
	private void loadRuntimes(File workDir) throws IllegalArgumentException, SecurityException, InstantiationException, IllegalAccessException, InvocationTargetException, NoSuchMethodException, ClassNotFoundException, IOException, SAXException, ParserConfigurationException  { 
		logger.info("loading runtimes");

		File runtimeConfig = new File(DirectoryLayout.getInstance().getConfDir(), "runtimes.xml");

		Document document = XmlUtil.parseReader(new FileReader(runtimeConfig));
		Element runtimesElement = (Element)document.getElementsByTagName("runtimes").item(0);


		for (Element runtimeElement: XmlUtil.getChildElements(runtimesElement, "runtime")) {
			String runtimeName = XmlUtil.getChildElement(runtimeElement, "name").getTextContent();
			logger.info("loading runtime " + runtimeName);
			if (runtimes.containsKey(runtimeName)) {
				logger.warn("runtime with the same name " + runtimeName + " already loaded, keeping the first one");
				continue;
			}
			
			
			boolean runtimeDisabled = runtimeElement.getAttribute("disabled").equals("true");
			if (runtimeDisabled) {
				logger.info("runtime " + runtimeName + " disabled by config");
			}
			
			// handler
			Element handlerElement = XmlUtil.getChildElement(runtimeElement, "handler");
			String handlerClassName = XmlUtil.getChildElement(handlerElement, "class").getTextContent();

			// parameters to handler
			HashMap<String, String> parameters = new HashMap<String, String>();

			// comp work dir
			parameters.put("workDir", workDir.toString());

			// parameters from config
			for (Element parameterElement: XmlUtil.getChildElements(handlerElement, "parameter")) {
				String paramName = XmlUtil.getChildElement(parameterElement, "name").getTextContent();
				String paramValue = XmlUtil.getChildElement(parameterElement, "value").getTextContent(); 
				parameters.put(paramName, paramValue);
			}

			// instantiate handler
			AnalysisHandler handler = (AnalysisHandler)Class.forName(handlerClassName).getConstructor(HashMap.class).newInstance(parameters);
			
			// disabled
			if (handler.isDisabled()) {
				runtimeDisabled = true;
				logger.info("runtime " + runtimeName + " disabled as handler is disabled");
			}

			// add to runtimes
			ToolRuntime runtime = new ToolRuntime(runtimeName, handler, runtimeDisabled); 
			this.runtimes.put(runtimeName, runtime);
		}
	}	

	// TODO: if we use modules we can probably remove this for now since
	//       we send information only about one module to the client.
	private void loadModules() throws IOException, SAXException, ParserConfigurationException { 
		logger.info("loading modules");
		
		String [] moduleFiles = new String[] { "microarray-module.xml" };
		for (String toolFileName: moduleFiles) {
			File moduleFile = new File(DirectoryLayout.getInstance().getConfDir(), toolFileName);
			if (moduleFile.exists()) {
				logger.info("loading from " + toolFileName);
				loadModule(moduleFile);
			}
		}
	}

	/**
	 * Parses a module file and loads all tools listed in it.
	 * 
	 * During parsing also construct description message that
	 * can be sent to the client.
	 * 
	 * @param toolFile
	 * @throws FileNotFoundException
	 * @throws SAXException
	 * @throws IOException
	 * @throws ParserConfigurationException
	 */
	private void loadModule(File toolFile) throws FileNotFoundException, SAXException,
	                                                     IOException, ParserConfigurationException {
		File toolConfig = toolFile;

		Document document = XmlUtil.parseReader(new FileReader(toolConfig));
		Element moduleElement = (Element)document.getElementsByTagName("module").item(0);
		
		// Construct a description message
		descriptionMessage = new DescriptionMessage(moduleElement.getAttribute("name"));
		
		// Stats
	    int totalCount = 0;
		int successfullyLoadedCount = 0;
	    int hiddenCount = 0;
	    int disabledCount = 0;
		
		// Read cetegories
		for (Element categoryElement: XmlUtil.getChildElements(moduleElement, "category")) {
		    
		    // Category's visibility
		    boolean categoryHidden = categoryElement.getAttribute("hidden").equals("true");
		    
		    if (!categoryHidden) {
		          
	            // Add to message
		        Category category = new Category(categoryElement.getAttribute("name"),
                                                 categoryElement.getAttribute("color"));		        
	            // Read tools
    		    for (Element toolElement: XmlUtil.getChildElements(categoryElement, "tool")) {
    		        totalCount++;
    		        
    		        String descriptionFilename = toolElement.getTextContent();
                    logger.debug("loading " + descriptionFilename);
                    
                    // Tool's visibility
                    boolean toolDisabled = toolElement.getAttribute("disabled").equals("true");
                    boolean toolHidden = toolElement.getAttribute("hidden").equals("true");
                    
                    // Runtime
                    String runtimeName = toolElement.getAttribute("runtime");
                    ToolRuntime runtime = runtimes.get(runtimeName);
                    if (runtime == null) {
                        logger.warn("loading " + descriptionFilename + " failed, could not find runtime " + runtimeName);
                        continue;
                    }
                    
                    // Description
                    AnalysisDescription description;
                    try {
                        description = runtime.getHandler().handle(descriptionFilename);
                    } catch (Exception e) {
                        logger.warn("loading " + descriptionFilename + " failed, could not create description", e);
                        continue;
                    }
                    if (descriptions.containsKey(description.getFullName())) {
                        logger.warn("loading " + descriptionFilename + " failed, description with the name " +
                                    description.getFullName() + " already exists");
                        continue;
                    }
                    descriptions.put(description.getFullName(), description);
                    
                    // Add to message
                    category.addTool(description.getName(), description.getSADL());
                    
                    successfullyLoadedCount++;
                    
                    // Supported and visible description lists
                    String disabledStatus = "";
                    if (!runtime.isDisabled() && !toolDisabled) {
                        supportedDescriptions.put(description.getFullName(), description);
                    } else {
                        disabledStatus = " DISABLED";
                        disabledCount++;
                    }
                    String hiddenStatus = "";
                    if (!toolHidden) {
                        visibleDescriptions.put(description.getFullName(), description);
                    } else {
                        hiddenStatus = " HIDDEN";
                        hiddenCount++;
                    }

                    logger.info("loaded " + description.getFullName().replace("\"", "") + " " +
                                description.getSourceResourceFullPath() + disabledStatus + hiddenStatus);
    		    }
    		    
                // Add to message after the tools have been read
                descriptionMessage.addCategory(category);
		    }
		}
		
        logger.info("loaded " + successfullyLoadedCount + "/" + totalCount +
                    " tools, " + disabledCount + " disabled, " + hiddenCount + " hidden");
	}

}