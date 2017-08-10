﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace principia {
namespace tools {

class StackTraceDecoder {
  const string dbh =
      @"\Program Files (x86)\Windows Kits\10\Debuggers\x64\dbh.exe";

  private static void Main(string[] args) {
    bool unity_crash;
    string commit = null;
    if (args.Length == 2) {
      unity_crash = false;
    } else if (args.Length == 3) {
      var match = Regex.Match(args[2],
                              "--unity-crash-at-commit=([0-9a-f]{40})");
      if (match.Success) {
        unity_crash = true;
        commit = match.Groups[1].ToString();
      } else {
        PrintUsage();
        return;
      }
    } else {
      PrintUsage();
      return;
    }
    string info_file_uri = args[0];
    string pdb_file = args[1];
    var web_client = new WebClient();
    var stream = new StreamReader(web_client.OpenRead(info_file_uri),
                                  Encoding.UTF8);
    if (!unity_crash) {
      var version_regex = new Regex(
          @"^I.*\] Principia version " + 
          @"([0-9]{10}-\w+)-[0-9]+-g([0-9a-f]{40}) built");
      Match version_match;
      do {
        version_match = version_regex.Match(stream.ReadLine());
      } while (!version_match.Success);
      string tag = version_match.Groups[1].ToString();
      commit = version_match.Groups[2].ToString();
    }
    var base_address_regex = new Regex(
        unity_crash ? @"GameData\\Principia\\principia.dll:principia.dll " +
                      @"\(([0-9A-F]+)\)"
                    : @"^I.*\] Base address is ([0-9A-F]+)$");
    Match base_address_match;
    do {
      base_address_match = base_address_regex.Match(stream.ReadLine());
    } while (!base_address_match.Success);
    string base_address_string = base_address_match.Groups[1].ToString();
    Int64 base_address = Convert.ToInt64(base_address_string, 16);
    Console.WriteLine("<!--- Using base address " +
                      Convert.ToString(base_address, 16) + " -->");
    var stack_regex = new Regex(
        unity_crash ? @"\(0x([0-9A-F]+)\) .*"
                    : @"@\s+[0-9A-F]+\s+.* \[0x([0-9A-F]+)(\+[0-9]+)?\]");
    Match stack_match;
    if (unity_crash) {
      Match stack_start_match;
      do {
        stack_start_match =
            Regex.Match(stream.ReadLine(),
                        @"========== OUTPUTING STACK TRACE ==================");
      } while (!stack_start_match.Success);
    }
    do {
      stack_match = stack_regex.Match(stream.ReadLine());
    } while (!stack_match.Success);
    var file_regex = new Regex(
        @"file\s+:\s+.*\\principia\\([a-z_]+)\\(\S+)");
    var line_regex = new Regex(@"line\s+:\s+([0-9]+)");
    for (;
         stack_match.Success;
         stack_match = stack_regex.Match(stream.ReadLine())) {
      Int64 address = Convert.ToInt64(stack_match.Groups[1].ToString(), 16);
      Int64 dbh_base_address = 0x1000000;
      string rebased_address =
          Convert.ToString(address - base_address + dbh_base_address, 16);
      var p = new Process();
      p.StartInfo.UseShellExecute = false;
      p.StartInfo.RedirectStandardOutput = true;
      p.StartInfo.FileName = dbh;
      p.StartInfo.Arguments =
          '"' + pdb_file + "\" laddr \"" + rebased_address + '"';
      p.Start();
      string output = p.StandardOutput.ReadToEnd();
      p.WaitForExit();
      Match file_match = file_regex.Match(output);
      if (file_match.Success) {
        string file = file_match.Groups[1].ToString() + '/' +
                      file_match.Groups[2].ToString();
        string line = line_regex.Match(output).Groups[1].ToString();
        string url = "https://github.com/mockingbirdnest/Principia/blob/" +
                     commit + '/' + file + "#L" + line;
        Console.WriteLine("[`" + file + ":" + line + "`](" + url + ")");
      } else {
        Console.WriteLine("<!--- Nothing for " + stack_match.Groups[0] +
                          " -->");
      }
    }
  }

  private static void PrintUsage() {
    Console.WriteLine("Usage: stacktrace_decoder " +
                      "<info_file_uri> <principia_pdb_file> " +
                      "[--unity-crash-at-commit=<sha1>]");
  }
}

}  // namespace tools
}  // namespace principia
