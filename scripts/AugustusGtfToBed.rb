#!/usr/bin/env ruby
require 'ostruct'

FIELDS = %w{chrom source type start stop score strand frame attributes}

def lines_to_bed(transcript_id, lines)
  lines = lines.sort_by{|l| l.start.to_i}
  out = []
  out << lines.first.chrom
  out << lines.map{|l| l.start.to_i}.min - 1
  out << lines.map{|l| l.stop.to_i}.max
  out << transcript_id
  out << '.'
  out << lines.first.strand
  out << lines.map{|l| l.start.to_i}.min - 1
  out << lines.map{|l| l.stop.to_i}.max
  out << "."
  out << lines.count
  out << lines.sort_by{|l| l.start.to_i}.map{|l| l.stop.to_i - l.start.to_i + 1}.join(',')
  start_position = lines.map{|l| l.start.to_i}.min
  out << lines.sort_by{|l| l.start.to_i}.map{|l| l.start.to_i - start_position}.join(',')
  out.join("\t")
end

ARGF
  .map{|line| line.chomp}
  .find_all{|line| line.split("\t")[2] == "CDS"}
  .map{|line| OpenStruct.new(Hash[FIELDS.zip(line.split("\t"))])}
  .group_by{|line| match = line.attributes.match(/transcript_id \"([^\"]+)\"/); match.nil? ? nil : match[1]}
  .map{|transcript_id, lines| lines_to_bed(transcript_id, lines)}
  .each do |e|
  puts e
end
