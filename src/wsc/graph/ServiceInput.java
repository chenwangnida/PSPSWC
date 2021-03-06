package wsc.graph;

public class ServiceInput {

	private String input;
	boolean isSatified;

	public ServiceInput(String input, boolean isSatified) {
		super();
		this.input = input;
		this.isSatified = isSatified;
	}

	public String getInput() {
		return input;
	}

	public void setInput(String input) {
		this.input = input;
	}

	public boolean isSatified() {
		return isSatified;
	}

	public void setSatified(boolean isSatified) {
		this.isSatified = isSatified;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((input == null) ? 0 : input.hashCode());
		result = prime * result + (isSatified ? 1231 : 1237);
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ServiceInput other = (ServiceInput) obj;
		if (input == null) {
			if (other.input != null)
				return false;
		} else if (!input.equals(other.input))
			return false;
		if (isSatified != other.isSatified)
			return false;
		return true;
	}

}
