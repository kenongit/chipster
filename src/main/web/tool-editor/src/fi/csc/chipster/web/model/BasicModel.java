package fi.csc.chipster.web.model;

import com.vaadin.ui.Component;
import com.vaadin.ui.GridLayout;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TextArea;
import com.vaadin.ui.TextField;

import fi.csc.microarray.description.SADLDescription.Name;

public abstract class BasicModel extends GridLayout{
	
	private static final long serialVersionUID = 7942793110778486624L;
	public static final String OPTIONAL = "optional";
	public static final String NOT_OPTIONAL = "not optional";
	public static final String SINGLE_FILE = "Single file";
	public static final String MULTI_FILE = "Multiple files";
	public static final String MULTI_FILE_TEXT = "{...}";
	
	protected Label lbName;
	protected Label lbId;
	protected Label lbDescription;
	
	
	protected TextField name;
	protected TextField id;
	protected TextArea description;
	protected TextField prefix;
	protected TextField postfix;
	
	protected HorizontalLayout layout;
	
//	protected GridLayout grid;
	
	public BasicModel() {
//		grid = new GridLayout();
		this.setColumns(2);
		this.setImmediate(true);
		this.setMargin(true);
		this.setSpacing(true);
		this.setColumnExpandRatio(0, 2);
//		grid.setColumnExpandRatio(1, 2);
//		grid.setWidth("100%");
//		grid.setHeight("100%");
		initElements();
		addRow(lbId, id);
		addRow(lbName, name);
//		grid.addComponent(name, 0, 0);
//		
//		grid.addComponent(id, 1,0);
//		grid.addComponent(description, 2, 0, 2, 1);
	}
	
	private void initElements() {
		
		lbName = new Label("Display name:");
		lbId = new Label();
		lbDescription = new Label("Description");
		
		name = new TextField();
		name.setWidth("200px");
		
		id = new TextField();
		id.setWidth("200px");
		
		description = new TextArea();
		description.setWidth("200px");
		
		layout = new HorizontalLayout();
		
		prefix = new TextField();
		postfix = new TextField();
		
		layout.addComponent(prefix);
		layout.addComponent(new Label(MULTI_FILE_TEXT));
		layout.addComponent(postfix);
	}
	
	protected void addRow(Component label, Component component) {
		this.addComponent(label);
		this.addComponent(component);
	}
	
	protected String getValue(String value) {
		return value == null || value.isEmpty() ? "" : value;
	}
	
	protected String getValueOrNull(String value) {
		return value == null || value.isEmpty() ? null : value;
	}
	
	public abstract GridLayout createUI();
	
	protected void getMultipleFilesUI() {
		Area area = this.getComponentArea(id);
		if(area == null)
			return;
		this.removeComponent(id);
		this.addComponent(layout, area.getColumn1(), area.getRow1());
		
	}
	
	protected void getSingleFileUI() {
		Area area = this.getComponentArea(layout);
		if(area == null)
			return;
		this.removeComponent(layout);
		this.addComponent(id, area.getColumn1(), area.getRow1());
	}
	
	protected Name getNameFromUI(String value) {
		if(value.equals(SINGLE_FILE)) {
			return Name.createName(id.getValue(), name.getValue());
		} else if(value.equals(MULTI_FILE)) {
			return Name.createNameSet(prefix.getValue(), postfix.getValue(), name.getValue());
		} else {
			return null;
		}
	}

}
