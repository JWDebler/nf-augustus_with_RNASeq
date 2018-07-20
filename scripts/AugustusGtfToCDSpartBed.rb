#!/usr/bin/env ruby
require 'ostruct'

FIELDS = %w{chrom source type start stop score strand frame attributes}

ARGF
  .map{|line| line.chomp}
  .find_all{|line| line.split("\t")[2] == "CDSpart"}
  .map{|line| OpenStruct.new(Hash[FIELDS.zip(line.split("\t"))])}
  .each do |l|
  out = []
  out << l.chrom
  out << l.start.to_i - 1
  out << l.stop
  out << l.source
  out << l.score
  out << l.strand
  puts out.join("\t")
end

