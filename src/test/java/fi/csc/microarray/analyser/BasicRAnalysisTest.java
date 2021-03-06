/*
 * Created on Jan 31, 2005
 *
 */
package fi.csc.microarray.analyser;

import java.io.IOException;

import org.testng.annotations.Test;

import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.exception.MicroarrayException;

public class BasicRAnalysisTest extends AnalysisTestBase {

	@Test(groups = { "smoke"})
	public void testRNoOperation() throws Exception {
//		Task job = executor.createTask("\"Test\"/\"No-op\"");
//		executeJob(job);
	}

	@Test(groups = { "smoke"})
	public void testREcho() throws Exception {
//		Task job = executor.createTask("\"Test\"/\"Echo\"");
//		String input = "Haloo!";
//		
//		// FIXME rewrite job.addParameter("echoParameter", input);
//		executeJob(job);
//		
//		System.out.println("Output: " + job.getScreenOutput());
//		Assert.assertTrue(job.getScreenOutput().contains("[1] \"" + input + "\""));
	}
	
	@Test(groups = { "smoke"})
	public void testFailing() throws Exception {
//		Task job = executor.createTask("\"Test\"/\"Fail\"");
//		executeJob(job, State.FAILED);
//		Assert.assertTrue(job.getScreenOutput() != null);
	}
	
	@Test(groups = { "smoke"}, threadPoolSize = 1, invocationCount = 500)
	public void testInputOutput() throws Exception {
		
//		
//		Task job = executor.createTask("\"Test\"/\"InputOutput\"");
//		DataBean input = manager.createDataBean("input.tsv", new ByteArrayInputStream("input output test".getBytes()));
//		job.addInput("input.tsv", input);
//		
//		executeJob(job);
//		
//
//		DataBean output = job.getOutput("output.tsv");
//		
//		ByteArrayOutputStream out = new ByteArrayOutputStream();
//		IO.copy(output.getContentByteStream(), out);
//		
//		Assert.assertEquals("1\n2\n3\n", out.toString());
//		
	}

	@Test(groups = { "smoke"})
	public void testUniqueInputOutput() throws Exception {
//		
//		Task job = executor.createTask("\"Test\"/\"SimpleInputOutput\"");
//
//		String inputString = Thread.currentThread().getName() + "_" + System.currentTimeMillis() + "\n";
//
//		
//		
//		DataBean input = manager.createDataBean("input.dat", new ByteArrayInputStream(inputString.getBytes()));
//		job.addInput("input.dat", input);
//		
//		executeJob(job, 360, TimeUnit.SECONDS);
//		
//		DataBean output = job.getOutput("output.dat");
//		ByteArrayOutputStream out = new ByteArrayOutputStream();
//		
//		IO.copy(output.getContentByteStream(), out);
//		String outputString = new String(out.toByteArray());
//
//		String expectedString = "INPUT WAS: " + inputString;
//		if (!expectedString.equals(outputString)) {
//			System.out.println("EX: " + expectedString + " OUT: " + outputString);
//		}
//		
//		Assert.assertEquals(expectedString, outputString);
	}

	@Test(groups = { "smoke"})
	public void fetchDescriptions() throws TaskException, InterruptedException, IOException {
//		Task job = executor.createTask("describe");
//		executeJob(job);
//		String descriptions = new String(job.getOutput("description").getContents());
//		Assert.assertTrue(descriptions.contains("ANALYSIS"));
	}
	
	@Test(groups = { "smoke"})
	public void fetchSourceCode() throws TaskException, InterruptedException, IOException {
//		Task job = executor.createTask("describe-operation");
//		// FIXME rewrite job.addParameter("operation-name", "\"Normalisation\"/\"Affymetrix\"");
//		executeJob(job);
//		String sourceCode = new String(job.getOutput("sourcecode").getContents());
//		Assert.assertTrue(sourceCode.contains("Affymetrix normalization"));
	}
	
	
	@Test
	public void testInputSleepOutput() throws TaskException, InterruptedException, MicroarrayException {
//		Task task = executor.createTask("\"Test\"/\"InputSleepOutput\"");
//		String[] inputNames = new String[] {"GSM11805.cel", "GSM11814.cel", "GSM11823.cel", "GSM11830.cel"}; 
//		
//		int i = 1;
//		for (String name: inputNames) {
//			task.addInput("input" + i + ".cel", manager.createDataBean("input" + i + ".cel", this.getClass().getResourceAsStream("/kidney4-affy/" + name)));
//			i++;
//		}
//		executeJob(task);
	}
	

	@Test(groups = { "smoke"})
	public void testCancelTask() throws TaskException, InterruptedException, MicroarrayException {
//		Task task = executor.createTask("\"Test\"/\"InputSleepOutput\"");
//		String[] inputNames = new String[] {"affy_example.cel"};
//		
//		int i = 1;
//		for (String name: inputNames) {
//			task.addInput("input" + i + ".cel", manager.createDataBean("input" + i + ".cel", this.getClass().getResourceAsStream("/" + name)));
//			i++;
//		}
//		
//		CountDownLatch latch = new CountDownLatch(1);
//		task.addTaskEventListener(new JobResultListener(latch));
//		executor.startExecuting(task);
//		
//		Thread.sleep(1000);
//		
//		// cancel the job
//		executor.kill(task);
//
//		assert (task.getState() == State.CANCELLED);
	}
	
	@Test(threadPoolSize = 5, invocationCount = 30, groups = { "stress" })
	public void testRNoOperationStress() throws Exception {
		testREcho();
	}
	
	@Test(threadPoolSize = 10, invocationCount = 200, groups = { "stress"})
	public void testUniqueInputOutputStress() throws Exception {
		testUniqueInputOutput();
	}	
	
	@Test(groups = {"smoke"} )
	public void testScreenOutput() throws Exception {
//		Task job = executor.createTask("\"Test\"/\"Fail\"");
//		executeJob(job, State.FAILED);
//		System.out.println("Output: ");
//		System.out.println(job.getScreenOutput());
//		
//		Assert.assertTrue(job.getScreenOutput().startsWith("\n> failme"));
//		Assert.assertTrue(job.getScreenOutput().contains("Error:"));
	}

	@Test(groups = {"stress"}, threadPoolSize = 10, invocationCount = 100 )
	public void testScreenOutputStress() throws Exception {
//		Task job = executor.createTask("\"Test\"/\"Fail\"");
//		executeJob(job, State.FAILED);
//		System.out.println("Output: ");
//		System.out.println(job.getScreenOutput());
//		
//		Assert.assertTrue(job.getScreenOutput().startsWith("\n> failme"));
//		Assert.assertTrue(job.getScreenOutput().contains("Error:"));
	}	

}
