package fi.csc.microarray.analyser;

import java.io.File;
import java.io.FileInputStream;
import java.util.LinkedList;

import org.testng.Assert;
import org.testng.annotations.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.module.chipster.ChipsterSADLParser.Validator;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.XmlUtil;

public class SADLDescriptionTest {

	public static void main(String[] args) throws Exception {
		SADLDescriptionTest test = new SADLDescriptionTest();
		test.testDescriptions();
	}
	
	@Test(groups = {"unit"} )
	public void testDescriptions() throws Exception {
				
		DirectoryLayout.initialiseUnitTestLayout();
		
		LinkedList<String> resources = new LinkedList<String>();
		
		// Iterate through all tools and collect their resource definitions (filename, classname etc.) 
		for (File file : new File("src/main/applications/wrapper/comp/conf").listFiles()) {
			if (file.getName().endsWith("-module.xml")) {
				Document module = XmlUtil.parseFile(file);
				NodeList tools = module.getDocumentElement().getElementsByTagName("tool");
				for (int i = 0; i < tools.getLength(); i++) {
					Element tool = (Element)tools.item(i);
					String resource = tool.getElementsByTagName("resource").item(0).getTextContent();
					
					if (resource.endsWith(".acd")) {
						// Refers to EMBOSS ACD, they are converted in the code level, so there is not much to test
						continue;
					}

					resources.add(resource);
				}
			}
		}
		
		// Load SADL from each resource
		for (String resource : resources) {
			try {
				String sadl = null;
				
				// Try to figure out what it was. We have this information in the module descriptions,
				// but for keeping the test code simple, we infer it from the resource name. Currently
				// that is enough, in future we might need to use the actual module loading facility
				// to parse the module files.
				if (resource.split("\\.").length > 2) {
					// Is a class name

					System.out.println("validating class " + resource);
					JavaAnalysisJobBase jobBase = (JavaAnalysisJobBase)Class.forName(resource).newInstance();
					sadl = jobBase.getSADL();
					
				} else { 
					// Is a file name
					
					// Determine which file it is
					File[] dirsContainingDescriptions = new File[] {
						new File("src/main/tools/bsh"),
						new File("src/main/tools/shell"),
						new File("src/main/tools/R-2.12"),
						new File("src/main/tools/R-2.10"),
						new File("src/main/tools/R-2.9"),
					};
					LinkedList<File> potentialFiles = new LinkedList<File>(); 

					for (File dir : dirsContainingDescriptions) {
						for (File file : dir.listFiles()) {
							potentialFiles.add(file);
						}
					}

					for (File file : potentialFiles) {
						if (file.getName().endsWith(resource)) {

							// Found the file, determine the type and process it
							if (resource.endsWith(".R")) {

								// Is an R script
								SADLTool.ParsedScript res = new SADLTool().parseScript(new FileInputStream(file), "#");
								sadl = res.SADL;

							} else if (resource.endsWith(".bsh")) {

								// Is a BeanShell script
								SADLTool.ParsedScript res = new SADLTool().parseScript(new FileInputStream(file), "//");
								sadl = res.SADL;

							} else {

								// IS a plain SADL file
								sadl = Files.fileToString(file);
							}

							System.out.println("validating file " + file.getCanonicalFile());
							
							break; // we are done with this file
						}
					}
					
				}

				// Finally, validate the description
				if (sadl != null) {

					new Validator().validate(resource, sadl);
					
				} else {
					throw new RuntimeException("don't know what to do with: " + resource);
				}
				
			} catch (Exception e) {
				e.printStackTrace();
				Assert.fail("when parsing " + resource + ": " + e.getMessage() + " (" + e.getClass().getSimpleName() + ")");
			}
		}

	}
}
