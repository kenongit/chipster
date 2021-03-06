package fi.csc.microarray.messaging.message;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.log4j.Logger;

/**
 * Adds ability to attach parameters into a regular ChipsterMessage.
 * Parameters are accessed by their order (FIFO). There
 * are two types of parameters: named parameters, used for
 * system level messaging, and anonymous (normal) parameters,
 * given by the user. The latter have tight access policy to
 * make sure that required safety precautions are in place.
 * Implementation of safety precautions is subclass specific,
 * but they should typically guard against injection attacks
 * and simple DOS attacks (too long String parameters etc.).
 * 
 * @author Aleksi Kallio
 */
public abstract class ParameterMessage extends ChipsterMessage {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(ParameterMessage.class);

	private static final String KEY_PARAMETER = "parameter";
	private static final String KEY_NAMED_PARAMETER_KEY ="named-parameter-key";
	private static final String KEY_NAMED_PARAMETER_VALUE ="named-parameter-value";
	
	public static final String PARAMETER_AS_ID = "as-id";
	public static final String PARAMETER_JOB_ID = "job-id";
	public static final String PARAMETER_USE_COMPRESSION = "use-compression";
	public static final String PARAMETER_DISK_SPACE = "disk-space";

	
	private List<String> parameters = new LinkedList<String>();
	private HashMap<String, String> namedParameters = new HashMap<String, String>();
	
	public ParameterMessage() {
	}
	
	public ParameterMessage(List<String> parameters) {
		if (parameters != null) {
			this.parameters.addAll(parameters);
		}
	}

	/**
	 * Gets parameters in the order they were inserted.
	 * Only available to subclasses. Subclasses must
	 * make sure that required safety precautions
	 * are in place.
	 */
	protected List<String> getParameters() {
		return parameters;
	}


	/**
	 * Adds a parameter.  
	 * @param parameter
	 */
	public void addParameter(String parameter) {
		parameters.add(parameter);
	}
	
	public void addNamedParameter(String key, String value) {
		namedParameters.put(key, value);
	}
	
	public String getNamedParameter(String key) {
		return namedParameters.get(key);
	}
	
	public Set<String> getParameterNames() {
		return namedParameters.keySet();
	}
	
	
	@Override
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
		
		// load parameters
		this.parameters = new ArrayList<String>();
		String input;
		for (int i = 0; from.itemExists(KEY_PARAMETER + i); i++) {
			input = from.getString(KEY_PARAMETER + i);
			logger.debug("parameter " + (KEY_PARAMETER + i) + " is " + new String(input));
			this.parameters.add(input);
		}

		// load named parameters
		String key;
		String value;
		for (int i = 0; from.itemExists(KEY_NAMED_PARAMETER_KEY + i); i++) {
			key = from.getString(KEY_NAMED_PARAMETER_KEY + i);
			value = from.getString(KEY_NAMED_PARAMETER_VALUE + i);
			logger.debug("Unmarshalled named parameter: " + key + " = " + value);
			this.namedParameters.put(key, value);
		}
		
	}
	
	@Override
	public void marshal(MapMessage to) throws JMSException {
		super.marshal(to);
		
		// add parameters
		int i = 0;
		for (String parameter : parameters) {
			logger.debug("populated map message with " + KEY_PARAMETER + i + ": " + new String(parameter));
			to.setString(KEY_PARAMETER + i, parameter);
			i++;
		}

		// add named parameters
		
		logger.debug("Marshalling named parameters: " + namedParameters.keySet());
		int j = 0;
		for (String key : namedParameters.keySet()) {
			to.setString(KEY_NAMED_PARAMETER_KEY + j, key);
			to.setString(KEY_NAMED_PARAMETER_VALUE + j, namedParameters.get(key));
			j++;
		}
	}
}
