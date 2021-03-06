package fi.csc.microarray.client.visualisation.methods;

import javax.swing.JComponent;
import javax.swing.JPanel;

import org.icepdf.ri.common.SwingController;
import org.icepdf.ri.common.SwingViewBuilder;
import org.icepdf.ri.common.views.DocumentViewControllerImpl;

import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

public class PDFViewer extends Visualisation {

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return bean.isContentTypeCompatitible("application/pdf");
		
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {

		// build a controller
		SwingController controller = new SwingController();

		// Build a SwingViewFactory configured with the controller
		SwingViewBuilder factory = new SwingViewBuilder(controller);

		// Use the factory to build a JPanel that is pre-configured
		//with a complete, active Viewer UI.
		JPanel viewerComponentPanel = factory.buildViewerPanel();

		// Open a PDF document to view
		controller.openDocument(data.getContentByteStream(), data.getName(), data.getContentUrl().toString());

		// Set view mode
		controller.setPageViewMode(DocumentViewControllerImpl.ONE_COLUMN_VIEW, true);

		
		return viewerComponentPanel;
	}
}
