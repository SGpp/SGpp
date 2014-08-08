#!/usr/bin/env ruby

MAIN_EXEC = './main'

def error str
  $stderr.puts str
  exit 1
end

ARGV.each { |file|

  # Parse configuration file
  
  error "Invalid file #{file}." if not File.exists? file

  config = {
    "mode" => "LEARNER_ONLINE_SGD",
    "experimentDir" => "./results/exp01",
    "mainARFFFilename" => "data/cod_rna.train.bcnorm.arff",
    "testARFFFilename" => "data/cod_rna.test.bcnorm.arff",
    "dim" => "3",
    "level" => "3",
    "refinementCondition" => "FIXED_NUMBER",
    "refinementType" => "SURPLUS",
    "numIterations" => "10",
    "refinementNumPoints" => "2",
    "numRuns" => "1",
    "batchSize" => "10",
    "regularizationLambda" => "0.001",
    "CGStepSizeGamma" => "0.0001",
    "errorType" => "ACCURACY"
  }

  File.read(file).each_line{ |l|
    l.strip!
    next if l == ''

    l = l.split(":").map{|x|x.strip}
    config[l[0]] = l[1]
  }

  # Run main
  
  args = []
  config.each_pair{|k,v| 
    args << "#{k} #{v}" 
  }

  log_path = config["experimentDir"] + "/log.stdout"
  File.delete(log_path) if File.exists? log_path

  cmd = "#{MAIN_EXEC} #{args.join(" ")}"
  IO.popen(cmd) {|io|
    while (line = io.gets)
      puts line
      open(log_path, 'a') {|f|
        f.puts line
      }
    end
  }
}
