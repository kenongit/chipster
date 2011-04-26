package fi.csc.microarray.client.session;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.transform.stream.StreamSource;


import org.apache.log4j.Logger;
import org.mortbay.io.WriterOutputStream;
import org.xml.sax.SAXException;

import fi.csc.microarray.client.NameID;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.session.schema.DataType;
import fi.csc.microarray.client.session.schema.FolderType;
import fi.csc.microarray.client.session.schema.InputType;
import fi.csc.microarray.client.session.schema.LinkType;
import fi.csc.microarray.client.session.schema.NameType;
import fi.csc.microarray.client.session.schema.OperationType;
import fi.csc.microarray.client.session.schema.ParameterType;
import fi.csc.microarray.client.session.schema.SessionType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataBean.StorageMethod;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.IOUtils;

public class SessionLoader {
	
	private File sessionFile;
	private SessionType sessionType;
	
	private LinkedHashMap<String, DataFolder> folders = new LinkedHashMap<String, DataFolder>();
	private HashMap<DataFolder, FolderType> folderTypes = new HashMap<DataFolder, FolderType>();

	private LinkedHashMap<String, DataBean> dataBeans = new LinkedHashMap<String, DataBean>();
	private HashMap<DataBean, DataType> dataTypes = new HashMap<DataBean, DataType>();

	private LinkedHashMap<String, OperationRecord> operationRecords = new LinkedHashMap<String, OperationRecord>();
	private HashMap<OperationRecord, OperationType> operationTypes = new HashMap<OperationRecord, OperationType>();

	
	private DataManager dataManager;
	
	private static final Logger logger = Logger.getLogger(SessionLoader.class);
	
	
	public SessionLoader(File sessionFile) throws MicroarrayException {
		if (!UserSession.isValidSessionFile(sessionFile)) {
			throw new MicroarrayException("Not a valid session file.");
		}
		this.sessionFile = sessionFile;
		this.dataManager = Session.getSession().getDataManager();
	}
	
	
	public void loadSession() throws ZipException, IOException, JAXBException, SAXException {

		// parse metadata to jaxb classes
		parseMetadata();

		// create the basic objects from the jaxb classes 
		createFolders();
		createDataBeans();
		createOperations();

		// create the links between the objects
		linkOperationsToOutputs();
		linkDataItemChildren(dataManager.getRootFolder());
		linkDataBeans();
		linkInputsToOperations();
		
	}

	private void parseMetadata() throws ZipException, IOException, JAXBException, SAXException {
		ZipFile zipFile = null;
		try {
			// get the session.xml zip entry
			zipFile = new ZipFile(sessionFile);
			InputStream metadataStream = zipFile.getInputStream(zipFile.getEntry(UserSession.SESSION_DATA_FILENAME));

			// validate
			//ClientSession.getSchema().newValidator().validate(new StreamSource(metadataStream));
			
			// parse the metadata xml to java objects using jaxb
			Unmarshaller unmarshaller = UserSession.getJAXBContext().createUnmarshaller();
			unmarshaller.setSchema(UserSession.getSchema());
			NonStoppingValidationEventHandler validationEventHandler = new NonStoppingValidationEventHandler();
			unmarshaller.setEventHandler(validationEventHandler);
			this.sessionType = unmarshaller.unmarshal(new StreamSource(metadataStream), SessionType.class).getValue();
			
			if (validationEventHandler.hasEvents()) {
				throw new JAXBException("Invalid session file:\n" + validationEventHandler.getValidationEventsAsString());
			}
		}
		finally {
			// try to close all input streams from the zip file
			if (zipFile != null) {
				zipFile.close();
			}
		}
	}

	private void createFolders() {
		for (FolderType folderType : sessionType.getFolder()) {
			String name = folderType.getName();
			String id = folderType.getId();
			
			// check for unique id
			if (getDataItem(id) != null) {
				logger.warn("duplicate folder id: " + id + " , ignoring folder: " + name);
				continue;
			}
			
			// create the folder
			DataFolder dataFolder;
			if (UserSession.ROOT_FOLDER_ID.equals(id)) {
				dataFolder = dataManager.getRootFolder();
			} else {
				dataFolder = dataManager.createFolder(name);
			}
			folders.put(id, dataFolder);
			folderTypes.put(dataFolder, folderType);
		}
	}

	private void createDataBeans() {
		for (DataType dataType : this.sessionType.getData()) {
			String name = dataType.getName();
			String id = dataType.getId();
			
			// check for unique id
			if (getDataItem(id) != null) {
				logger.warn("duplicate data bean id: " + id + " , ignoring data bean: " + name);
				continue;
			}
			
			// create the data bean
			String storageMethodString = dataType.getStorageType();
			StorageMethod storageMethod = DataBean.StorageMethod.valueOf(storageMethodString);
			String urlString = dataType.getUrl();
			URL url = null;
			try {
				url = new URL(urlString);
			} catch (MalformedURLException e) {
				logger.warn("could not parse url: "  + urlString + " for data bean: " + name);
				continue;
			}
			
			
			DataBean dataBean = null;
			switch (storageMethod) {
			case LOCAL_SESSION:
				try {
					dataBean = dataManager.createDataBeanFromZip(name, url);
				} catch (MicroarrayException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				break;
			case LOCAL_USER:
				try {
					dataBean = dataManager.createDataBean(name, url);
				} catch (MicroarrayException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				break;
			default:
				logger.warn("unexpected storage method " + storageMethod.name() + " for data bean: " + name);	
				continue;
			}

			// cache url
			String cacheURLString = dataType.getCacheUrl();
			if (cacheURLString != null) {
				URL cacheURL = null;
				try {
					cacheURL = new URL(cacheURLString);
				} catch (MalformedURLException e1) {
					logger.warn("could not parse cache url: "  + cacheURLString + " for data bean: " + name);
				}
				dataBean.setCacheUrl(cacheURL);
			}

			// notes
			dataBean.setNotes(dataType.getNotes());
			//			dataBean.setCreationDate(date);
			
			dataBean.setContentType(Session.getSession().getDataManager().guessContentType(dataBean.getName()));
			
			dataBeans.put(id, dataBean);
			dataTypes.put(dataBean, dataType);
		}
	}

	
	private void createOperations() {
		for (OperationType operationType : sessionType.getOperation()) {
			String operationSessionId = operationType.getId();

			// check for unique id
			if (operationRecords.containsKey(operationSessionId)) {
				logger.warn("duplicate operation id: " + operationSessionId);
				continue;
			}

			OperationRecord operationRecord = new OperationRecord();

			// name id
			operationRecord.setNameID(createNameID(operationType.getName()));
			
			// category
			operationRecord.setCategoryName(operationType.getCategory());
			String colorString = operationType.getCategoryColor();
			if (colorString != null) {
				operationRecord.setCategoryColor(Color.decode(colorString));
			}
			
			// parameters
			for (ParameterType parameterType : operationType.getParameter()) {
				operationRecord.addParameter(parameterType.getName().getId(), parameterType.getName().getDisplayName(), parameterType.getName().getDescription(), parameterType.getValue());
			}

			// source code
			String sourceCodeFileName = operationType.getSourceCodeFile();
			if (sourceCodeFileName != null && !sourceCodeFileName.isEmpty()) {
				String sourceCode = null;
				try {
					sourceCode = getSourceCode(sourceCodeFileName);
				} catch (Exception e) {
					logger.warn("could not load source code from " + sourceCodeFileName);
				}

				// might be null, is ok
				operationRecord.setSourceCode(sourceCode);
			}
			
			operationRecords.put(operationSessionId, operationRecord);
			operationTypes.put(operationRecord, operationType);
		}
	}

	
	private void linkDataItemChildren(DataFolder parent) {
		for (String childId : folderTypes.get(parent).getChild()) {
			
			// check that the referenced data item exists
			DataItem child = getDataItem(childId);
			if (child == null) {
				logger.warn("child with id: " + childId + " does not exist");
				continue;
			}

			// add as a child
			parent.addChild(child);
			
			// recursively go inside folders
			if (child instanceof DataFolder) {
				linkDataItemChildren((DataFolder) child);
			}
		}
	}
	
	/**
	 * Link OperationRecords and DataBeans by adding real input DataBean references
	 * to OperationRecords.
	 * 
	 */
	private void linkInputsToOperations() {
		
		for (OperationRecord operationRecord : operationRecords.values()) {
			
			// get data bean ids from session data
			for (InputType inputType : operationTypes.get(operationRecord).getInput()) {
				
				String inputID = inputType.getData();
				if (inputID == null) {
					continue;
				}
				
				// find the data bean
				DataBean inputBean = dataBeans.get(inputID);
				if (inputBean == null) {
					continue;
				}
				
				// skip phenodata, it is bound automatically
				if (inputBean.queryFeatures("/phenodata/").exists()) {
					continue; 
				}

				// add the reference to the operation record
				operationRecord.addInput(createNameID(inputType.getName()), inputBean);
			}
		}
	}

	
	/**
	 * Add links form DataBeans to the OperationRecord which created the DataBean.
	 * 
	 * If OperationRecord is not found, use unknown OperationRecord.
	 * 
	 */
	private void linkOperationsToOutputs() {
		
		for (DataBean dataBean : dataBeans.values()) {
			String operationId = dataTypes.get(dataBean).getResultOf();
			OperationRecord operationRecord = null;
			if (operationId != null) {
				operationRecord = operationRecords.get(operationId);
			}

			// if operation record is not found use dummy
			if (operationRecord == null) {
				operationRecord = OperationRecord.getUnkownOperationRecord();
			}
			
			dataBean.setOperationRecord(operationRecord);
		}
	}

	/**
	 * 
	 */
	private void linkDataBeans() {
		for (DataBean dataBean : dataBeans.values()) {
			for (LinkType linkType : dataTypes.get(dataBean).getLink()) {
				// if something goes wrong for this link, continue with others
				try {
					String targetID = linkType.getTarget();
					if (targetID != null) {
						DataBean targetBean = dataBeans.get(targetID);
						if (targetBean != null) {
							dataBean.addLink(Link.valueOf(linkType.getType()), targetBean);
						}
					}
					
				} catch (Exception e) {
					logger.warn("could not add link", e);
					continue;
				}
			}
			
		}
	}
	
	/**
	 * 
	 * @param id
	 * @return null if no DataItem for the id is found
	 */
	private DataItem getDataItem(String id) {
		DataItem dataItem = folders.get(id);
		if (dataItem != null) {
			return dataItem;
		} else { 
			return dataBeans.get(id);
		}
	}

	private NameID createNameID(NameType name) {
		return new NameID(name.getId(), name.getDisplayName(), name.getDescription());
	}

	private String getSourceCode(String sourceCodeFileName) throws ZipException, IOException {
		ZipFile zipFile = null;
		InputStream sourceCodeInputStream = null;
		StringWriter stringWriter = null;
		try {
			zipFile = new ZipFile(sessionFile);
			sourceCodeInputStream = zipFile.getInputStream(zipFile.getEntry(sourceCodeFileName));

			stringWriter = new StringWriter();
			IOUtils.copy(sourceCodeInputStream, new WriterOutputStream(stringWriter));
			stringWriter.flush();
		}
		finally {
			IOUtils.closeIfPossible(sourceCodeInputStream);
			IOUtils.closeIfPossible(stringWriter);
			if (zipFile != null) {
				zipFile.close();
			}
		}

		return stringWriter.toString();
	}

}
