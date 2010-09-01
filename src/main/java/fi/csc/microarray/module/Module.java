package fi.csc.microarray.module;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.List;

import javax.swing.JMenu;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.QuickLinkPanel;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;

/**
 * Client side module. Encapsulates all application area specific logic, e.g., DNA microarray module
 * encapsulates all DNA microarray specific functionality.
 * 
 * @author Aleksi Kallio
 *
 */
public interface Module {

	/**
	 * Plugs features of this module to given data manager.
	 * 
	 * @param manager data manager to plug into
	 */
	public void plugFeatures(DataManager manager);
	
	/**
	 * Plugs modifiers of this module to given data manager.
	 * 
	 * @param manager data manager to plug into
	 */
	public void plugModifiers(DataManager manager);
	
	/**
	 * Plugs content types of this module to given data manager.
	 * 
	 * @param manager data manager to plug into
	 */
	public void plugContentTypes(DataManager manager);
	
	/**
	 * Plugs type tags of this module to given data manager.
	 * 
	 * @param manager data manager to plug into
	 */
	public void plugTypeTags(DataManager manager);

	/**
	 * Returns the name of the server module associated to this 
	 * this client side module, or null if not available.
	 *  
	 * @return server module name or null
	 */
	public String getServerModuleName();
	
	/**
	 * Adds import menu items to menu.
	 * 
	 * @param importMenu menu to add links to
	 */
	public void addImportMenuItems(JMenu importMenu);
	
	/**
	 * Adds import links link list.
	 * 
	 * @param quickLinkPanel quick link panel used when creating the links
	 * @param importLinks link list to add to
	 */
	public void addImportLinks(QuickLinkPanel quickLinkPanel, List<JXHyperlink> importLinks);
	
	/**
	 * Is import tool supported when this module is used? 
	 * 
	 * @return true iff import tool is supported by this module
	 */
	public boolean isImportToolSupported();
	
	/**
	 * Checks data bean for workflow compatibility.
	 * 
	 * @param data data bean to check
	 * @return true iff data bean supports workflows under, when this module is used
	 */
	public boolean isWorkflowCompatible(DataBean data);
	
	/**
	 * Returns visualisation methods of this module.
	 * 
	 * @return array of visualisation methods
	 */
	public VisualisationMethod[] getVisualisationMethods();
	
	/**
	 * Return URL to example session file or null, if not available.
	 * @return url or null
	 * @throws MalformedURLException
	 */
	public URL getExampleSessionUrl() throws MalformedURLException;
}
