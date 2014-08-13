#!/usr/bin/ruby
#
# Usage:
# generate-exp.rb DATASET

require 'shellwords'

if ARGV.empty?
  ARGV=["breast-cancer", "cod-rna", "german", "ijcnn1", "magic04"]
end

ARGV.each {|name|
  puts "Process: #{name}"

  common_path="generate-exp/#{name}/common"
  var_path="generate-exp/#{name}/var"

  if (not File.exists? common_path) or (not File.exists? var_path)
    $stderr.puts "One of the configuration files doesn't exist."
    exit 1
  end

  common_content = File.read(common_path).strip

  File.read(var_path).split("\n\n").compact.each_with_index {|x, i|
    dir = "exp-#{name}-#{sprintf "%02d", (i+1)}"
    file = dir + "/config"
    var = x.gsub("{DIRNAME}", dir).strip
    all = common_content.gsub("{DIRNAME}", dir) + "\n" + var

    # puts "Create dir: #{dir}"
    `mkdir -p #{dir.shellescape}`
    puts "Write: #{file}"
    File.write(file, all)
  }
}
