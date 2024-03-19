version 1.0

workflow MHap_Analysis {
	call MHap_Analysis_tool
}

task MHap_Analysis_tool {
  input {
    String name
  }

  command {
    echo "Hello this is my ${name}"
  }
  
  output {
  
    File MHap_Analysis_output = stdout()
  
  }

}
