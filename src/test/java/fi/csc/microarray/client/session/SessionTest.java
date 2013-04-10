package fi.csc.microarray.client.session;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.LinkedList;
import java.util.Random;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.client.RemoteServiceAccessor;
import fi.csc.microarray.client.ServiceAccessor;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.messaging.MessagingTestBase;
import fi.csc.microarray.module.ModuleManager;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class SessionTest extends MessagingTestBase {

	
	public SessionTest() {
		super("demo", "ddemo");
	}
	
	@Test(groups = {"unit"} )
	public void testStorageSessions() throws Exception {
		
		// set up modules
		ModuleManager moduleManager = new ModuleManager("fi.csc.microarray.module.chipster.MicroarrayModule");
		Session.getSession().setModuleManager(moduleManager);
		
		// set up system
		DataManager manager = new DataManager();
		moduleManager.plugAll(manager, null);
		LinkedList<ToolModule> toolModules = new LinkedList<ToolModule>();
		ServiceAccessor serviceAccessor = new RemoteServiceAccessor();
		serviceAccessor.initialise(manager, authenticationListener);
		serviceAccessor.fetchDescriptions(new MicroarrayModule());
		Session.getSession().setServiceAccessor(serviceAccessor);
		toolModules.addAll(serviceAccessor.getModules());
		FileBrokerClient fileBrokerClient = serviceAccessor.getFileBrokerClient();
		
		// create data
		DataBean data = manager.createDataBean("test");
		data.setParent(manager.getRootFolder());
		int dbCountOrig = manager.databeans().size();
		
		// save storage session
		String sessionName = "unit test session " + new Random().nextInt(10000);
		manager.saveStorageSession(sessionName);
				
		// delete all data
		manager.deleteAllDataItems();
		
		// load session
		String[][] sessions = fileBrokerClient.listRemoteSessions();
		URL sessionURL = findSession(sessionName, sessions);
		manager.loadStorageSession(sessionURL);

		// assert loaded data is ok
		Assert.assertEquals(manager.databeans().size(), dbCountOrig);
		
		// remove remote session
		fileBrokerClient.removeRemoteSession(sessionURL);

		// assert removal
		String[][] sessions2 = fileBrokerClient.listRemoteSessions();
		Assert.assertNull(findSession(sessionName, sessions2));

	}

	private URL findSession(String sessionName,	String[][] sessions) throws MalformedURLException {
		URL sessionURL = null;
		for (int i = 0; i < sessions[0].length; i++) {
			if (sessionName.equals(sessions[0][i])) {
				sessionURL = new URL(sessions[1][i]);
				break;
			}
		}
		return sessionURL;
	}

}
