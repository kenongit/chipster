package fi.csc.microarray.client.dialog;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JCheckBox;
import javax.swing.JPanel;

/**
 * Accessory component for the JFileChooser to allow direct input of dataset
 * folder.
 * 
 * @author Petri Klemelä
 * 
 */
public class SnapshotAccessory extends JPanel implements ActionListener {
	
	private JCheckBox clearCheckBox;

	public SnapshotAccessory() {

		this.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weightx = 1;
		
		c.gridx = 0;
		c.gridy = 0;
		c.insets.set(10, 5, 0, 0);
		c.weighty = 1;
		clearCheckBox = new JCheckBox("<html>Add to current<br>session</html>");
		clearCheckBox.setSelected(false);
		this.add(clearCheckBox, c);	
	}

	public boolean clearSession() {
		return !clearCheckBox.isSelected();
	}

	public void setDefaults() {
		clearCheckBox.setSelected(false);
	}

	public void actionPerformed(ActionEvent e) {
	}
}