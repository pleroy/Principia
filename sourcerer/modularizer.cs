using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using static principia.sourcerer.Analyser;
using static principia.sourcerer.Filenames;
using static principia.sourcerer.Parser;
using static principia.sourcerer.Rewriter;

namespace principia {
namespace sourcerer {

class Modularizer {
  public static void ProcessProject(DirectoryInfo project,
                                    HashSet<string> excluded,
                                    HashSet<string> extra_headers,
                                    bool dry_run) {
        // Parse all the files in our project.
        var file_info_to_file = new Dictionary<FileInfo, Parser.File>();
    var file_name_to_file_info = new Dictionary<string, FileInfo>();
    FileInfo[] hpp_files = project.GetFiles("*.hpp");
    FileInfo[] cpp_files = project.GetFiles("*.cpp");
    FileInfo[] body_hpp_files = Array.FindAll(hpp_files, IsBodyHpp);
    FileInfo[] all_body_files = [.. body_hpp_files.Union(cpp_files)];
    FileInfo[] all_files = [.. hpp_files.Union(cpp_files)];
    foreach (FileInfo input_file in all_files) {
      if (excluded.Contains(input_file.Name)) {
        continue;
      }
      Parser.File parser_file =
          Parser.ParseFile(input_file, IsBody(input_file, extra_headers));
      file_info_to_file.Add(input_file, parser_file);
      file_name_to_file_info.Add(Path.GetFullPath(input_file.Name,
                                                  input_file.DirectoryName!),
                                 input_file);
    }

    // Map the bodies to their headers.
    var corresponding_header = new Dictionary<FileInfo, FileInfo>();
    foreach (FileInfo input_file in all_body_files) {
      string corresponding_header_file_name =
          Filenames.CorrespondingHeader(input_file);
      // Benchmarks and tests don't have a corresponding header.
      if (file_name_to_file_info.TryGetValue(corresponding_header_file_name,
                                             out FileInfo?
                                                     corresponding_header_file)) {
        if (corresponding_header_file != input_file) {
          corresponding_header.Add(input_file, corresponding_header_file);
        }
      }
    }

    // Fix our body files.
    var needs_import_std = new HashSet<FileInfo>();
    var body_principia_includes =
        new Dictionary<FileInfo, SortedSet<Include>>();
    var body_third_party_includes =
        new Dictionary<FileInfo, SortedSet<Include>>();
    foreach (FileInfo body_file in all_body_files) {
      if (excluded.Contains(body_file.Name) ||
          !IsBody(body_file, extra_headers)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[body_file];

      // Only process includes up to the first outermost namespace.  We don't
      // want to touch the final include (for the body) or includes of generated
      // code.
      var outermost_namespaces = FindOutermostNamespaces(parser_file);
      Namespace? first_namespace =
          outermost_namespaces.Count > 0 ? outermost_namespaces[0] : null;
      var includes = FindIncludes(parser_file, stop: first_namespace);

      // No link to our own header for .cpp files, the link is provided by the
      // module declaration.
      Node? link = null;
      if (!IsBodyCpp(body_file)) {
        link = FixOwnHeaderLink(parser_file, includes);
      }

      // Find all the includes and categorize them.
      SplitIncludes(includes,
                    excluded,
                    out SortedSet<Include> principia_includes,
                    out SortedSet<Include> system_includes,
                    out SortedSet<Include> third_party_includes,
                    out SortedSet<Include> conditional_includes);

      if (IsBodyHpp(body_file)) {
        FixPragmaOnce(parser_file, needs_export_module: false);
        // If this is a body that gets included into its header, its imports
        // will go in the header.
        FileInfo hpp_file = corresponding_header[body_file];
        if (!body_principia_includes.TryAdd(hpp_file, principia_includes)) {
          body_principia_includes[hpp_file].UnionWith(principia_includes);
        }
        if (system_includes.Count > 0) {
          needs_import_std.Add(hpp_file);
        }
        if (!body_third_party_includes.TryAdd(hpp_file, third_party_includes)) {
          body_third_party_includes[hpp_file].UnionWith(third_party_includes);
        }
      } else {
        Node? module;
        if (IsBodyCpp(body_file)) {
          module = FixModule(parser_file, after: link);
        } else {
          module = link;
        }
        Node? last_std_import;
        if (system_includes.Count > 0) {
          last_std_import = FixImportStd(parser_file, after: module);
        } else {
          last_std_import = module;
        }
        // If this is a test or benchmark, its imports go right here.
        Node? last_third_party_import =
            FixImportThirdPartyHeaders(parser_file,
                                       after: last_std_import,
                                       third_party_includes);
        FixImportPrincipiaHeaders(parser_file,
                                  after: last_third_party_import,
                                  principia_includes,
                                  skip_own_header: IsBodyCpp(body_file));
        FixImportConditionalHeaders(parser_file,
                                    conditional_includes);
      }

      // In all cases, the includes must go away.
      FixIncludes(parser_file, principia_includes);
      FixIncludes(parser_file, system_includes);
      FixIncludes(parser_file, third_party_includes);

      RewriteFile(body_file, parser_file, dry_run);
    }

    // Fix our non-body header files.  This must come last because we need the
    // information that we collected from the body files.
    foreach (FileInfo hpp_file in hpp_files) {
      if (excluded.Contains(hpp_file.Name) || IsBody(hpp_file, extra_headers)) {
        continue;
      }
      Parser.File parser_file = file_info_to_file[hpp_file];

      // Only process includes up to the first outermost namespace.  We don't
      // want to touch the final include (for the body) or includes of generated
      // code.
      var outermost_namespaces = FindOutermostNamespaces(parser_file);
      Namespace? first_namespace =
          outermost_namespaces.Count > 0 ? outermost_namespaces[0] : null;

      Node? last_export_module =
          FixPragmaOnce(parser_file, needs_export_module: true);

      SplitIncludes(FindIncludes(parser_file, stop: first_namespace),
                    excluded,
                    out SortedSet<Include> principia_includes,
                    out SortedSet<Include> system_includes,
                    out SortedSet<Include> third_party_includes,
                    out SortedSet<Include> conditional_includes);
      var all_principia_includes =
          new SortedSet<Include>(principia_includes, new IncludeComparer());
      all_principia_includes.UnionWith(
          body_principia_includes.GetValueOrDefault(
              hpp_file,
              new SortedSet<Include>(new IncludeComparer())));
      // Don't create an import for our own header, 'cause that would be
      // circular.
      all_principia_includes.RemoveWhere(inc => inc.is_own_header);
      var all_third_party_includes =
          new SortedSet<Include>(third_party_includes, new IncludeComparer());
      all_third_party_includes.UnionWith(
          body_third_party_includes.GetValueOrDefault(
              hpp_file,
              new SortedSet<Include>(new IncludeComparer())));

      // Add the imports collected from this file and its body.
      Node? last_std_import;
      if (system_includes.Count > 0 || needs_import_std.Contains(hpp_file)) {
        last_std_import = FixImportStd(parser_file, last_export_module);
      } else {
        last_std_import = last_export_module;
      }
      Node? last_third_party_import =
          FixImportThirdPartyHeaders(parser_file,
                                     after: last_std_import,
                                     all_third_party_includes);
      FixImportPrincipiaHeaders(parser_file,
                                after: last_third_party_import,
                                all_principia_includes,
                                skip_own_header: false);
      FixImportConditionalHeaders(parser_file,
                                  conditional_includes);

      // Drop the includes.
      FixIncludes(parser_file, principia_includes);
      FixIncludes(parser_file, system_includes);
      FixIncludes(parser_file, third_party_includes);

      // Export the non-internal using declarations.
      FixExportedUsingDeclarations(parser_file);

      // Export everything from the internal namespaces except the using
      // directives.
      FixInternalNamespacesExportSequence(parser_file);

      RewriteFile(hpp_file, parser_file, dry_run, move_to_ixx: true);
    }

    // Fix the project files to reference the .ixx files.
    FixProjectFile(
        new FileInfo(project.FullName + "\\" + project.Name + ".vcxproj"),
        extra_headers,
        dry_run);
    FixProjectFile(
        new FileInfo(
            project.FullName + "\\" + project.Name + ".vcxproj.filters"),
        extra_headers,
        dry_run);
  }

  // Ensures that `after` is followed by a blank line, unless it is `null` or
  // the last child of its parent.  Returns the blank line, or `null` or the
  // last child.
  private static Node? EnsureBlankLineAfter(Node? after) {
    if (after != null) {
      if (after.is_blank_line) {
        // If `after` is already a blank line we have nothing to do.
        return after;
      }
      Node parent = after.parent!;
      if (after.position_in_parent + 1 < parent.children.Count) {
        after = parent.children[after.position_in_parent + 1];
        if (!after.is_blank_line) {
          // If `after` is not followed by a blank line, insert one and move
          // `after` to that blank line.
          var blank = Text.Blank();
          parent.InsertChild(after.position_in_parent, blank);
          after = blank;
        }
      }
    }
    return after;
  }

  // Same as above, but for `before` and the preceding blank line.  Doesn't
  // return anything because we don't use it for sequencing the insertions.
  private static void EnsureBlankLineBefore(Node? before) {
    if (before != null) {
      if (before.is_blank_line) {
        // If `before` is already a blank line we have nothing to do.
        return;
      }
      Node parent = before.parent!;
      if (before.position_in_parent == 0 ||
          (before.position_in_parent - 1 >= 0 &&
           !parent.children[before.position_in_parent - 1].is_blank_line)) {
        // If `before` is not preceded by a blank line, insert one.
        var blank = Text.Blank();
        parent.InsertChild(before.position_in_parent, blank);
      }
    }
  }

  private static void FixExportedUsingDeclarations(Parser.File file) {
    var using_declarations = FindUsingDeclarations(file, internal_only: false);
    foreach (UsingDeclaration using_declaration in using_declarations) {
      Node parent = using_declaration.parent!;
      int position_in_parent = using_declaration.position_in_parent;
      var export = new Export(parent: null);
      parent.RemoveChild(using_declaration);
      export.AddChild(using_declaration);
      parent.InsertChild(position_in_parent, export);
    }
  }

  // Change conditional includes into imports, without moving them around.  They
  // may need fixing by hand afterwards.
  private static void FixImportConditionalHeaders(Parser.File file,
                                                  SortedSet<Include>
                                                      conditional_includes) {
    foreach (Include inc in conditional_includes) {
      if (inc.raw_path != null) {
        int position_in_parent = inc.position_in_parent;
        file.RemoveChild(inc);
        var import = new Import(parent: null, name: inc.raw_path);
        import.comment = inc.comment;
        file.InsertChild(position_in_parent, import);
      }
    }
  }

  // Inserts imports of the Principia headers in `principia_includes` after
  // `after` in `file`.  Returns the last inserted node.
  private static Node? FixImportPrincipiaHeaders(Parser.File file,
                                                 Node? after,
                                                 SortedSet<Include>
                                                     principia_includes,
                                                 bool skip_own_header) {
    int position_in_parent = (after?.position_in_parent + 1) ?? 0;
    Node? last = after;
    foreach (Include inc in principia_includes) {
      if (skip_own_header && inc.is_own_header) {
        continue;
      }
      var import =
          new Import(parent: null,
                     name: "principia." + string.Join('.', inc.path));
      import.comment = inc.comment;
      file.InsertChild(position_in_parent, import);
      ++position_in_parent;
      last = import;
    }
    return EnsureBlankLineAfter(last);
  }

  // Inserts an `import std` after `after` in `file` and returns it.
  private static Node? FixImportStd(Parser.File file,
                                    Node? after) {
    int position_in_parent = (after?.position_in_parent + 1) ?? 0;
    Node? last = null;
    var import_std = new Import(parent: null, "std");
    file.InsertChild(position_in_parent, import_std);
    last = import_std;

    // Surround the import with blank lines if needed.
    bool needs_blank_line_before =
        position_in_parent > 0 &&
        !file.children[position_in_parent - 1].is_blank_line;
    bool needs_blank_line_after =
        position_in_parent + 1 < file.children.Count &&
        !file.children[position_in_parent + 1].is_blank_line;

    if (needs_blank_line_before) {
      var blank = Text.Blank();
      file.InsertChild(import_std.position_in_parent, blank);
      last = blank;
    }
    if (needs_blank_line_after) {
      var blank = Text.Blank();
      file.InsertChild(import_std.position_in_parent + 1, blank);
      last = blank;
    }
    return last;
  }

  // Inserts imports of the third-party headers in `third_party_includes`
  // after `after` in `file`.  Returns the last inserted node.
  private static Node? FixImportThirdPartyHeaders(Parser.File file,
                                                  Node? after,
                                                  SortedSet<Include>
                                                      third_party_includes) {
    int position_in_parent = (after?.position_in_parent + 1) ?? 0;
    Node? last = after;
    foreach (Include inc in third_party_includes) {
      var import =
          new Import(parent: null,
                     name: '"' +
                           string.Join('/', inc.path) +
                           "."  +
                           inc.extension +
                           '"');
      import.comment = inc.comment;
      file.InsertChild(position_in_parent, import);
      ++position_in_parent;
      last = import;
    }
    return EnsureBlankLineAfter(last);
  }

  // Removes include directives for the given headers from the file.
  private static void
      FixIncludes(Parser.File file, SortedSet<Include> headers) {
    foreach (Include inc in headers) {
      int position_in_parent = inc.position_in_parent;
      file.RemoveChild(inc);
      // Remove blank lines that follow the include.
      while (position_in_parent < file.children.Count &&
             file.children[position_in_parent].is_blank_line) {
        file.RemoveChild(file.children[position_in_parent]);
      }
    }
  }

  private static void FixInternalNamespacesExportSequence(Parser.File file) {
    var internal_using_directives =
        FindUsingDirectives(file, internal_only: true);
    var es = new ExportSequence(parent: null);
    int position_in_parent = -1;
    Namespace? ns = null;
    if (internal_using_directives.Count == 0) {
      var innermost_namespaces = FindInnermostNamespaces(file);
      if (innermost_namespaces.Count == 0) {
        // No using directives and no namespaces, this is a weird file, leave it
        // alone.
        return;
      } else {
        // Insert the export sequence at the beginning of the first innermost
        // namespace.
        ns = innermost_namespaces[0];
        position_in_parent = 0;
      }
    } else {
      // Insert the export sequence after the last internal using directive.
      Node? after = internal_using_directives[^1];
      after = EnsureBlankLineAfter(after);
      ns = after!.parent as Namespace;
      position_in_parent = after.position_in_parent + 1;
    }
    ns!.InsertChild(position_in_parent, es);
    ++position_in_parent;
    while (position_in_parent < ns.children.Count) {
      Node node = ns.RemoveChild(position_in_parent);
      es.AddChild(node);
    }
    if (es.children.Count == 0) {
      // Nothing in the export sequence, remove it.
      ns.RemoveChild(es);
    } else {
      if (!es.children[0].is_blank_line) {
        // Start the sequence with a blank line.
        EnsureBlankLineBefore(es.children[0]);
      }
      if (!es.children[^1].is_blank_line) {
        // End the sequence with a blank line.
        EnsureBlankLineAfter(es.children[^1]);
      }
    }
  }

  // Adds a link to our own header if it is included in `file`, as the first
  // line of the file.  Returns it.
  private static Node? FixOwnHeaderLink(Parser.File file,
                                        List<Include> includes) {
    FileInfo file_info = file.file_info;
    foreach (Include inc in includes) {
      if (inc.is_own_header) {
        string usage =
            IsBenchmark(file_info) ? "Benchmark for" :
            IsTest(file_info) ? "Test for" :
            (IsBodyHpp(file_info) || IsBodyCpp(file_info)) ? "Body of" : "";
        var link =
            new Text(
                text: "// " +
                      usage +
                      " link:" +
                      string.Join('/', inc.path) +
                      ".ixx",
                parent: null);
        file.InsertChild(0, link);
        return EnsureBlankLineAfter(link);
      }
    }
    return null;
  }

  private static Node? FixModule(Parser.File file, Node? after) {
    int position_in_parent = (after?.position_in_parent + 1) ?? 0;
    var module = new Module(parent: null, name: file.file_module_full_name);
    file.InsertChild(position_in_parent, module);
    return EnsureBlankLineAfter(module);
  }

  // Removes `#pragma once` directives from the file and insert an
  // `export module` directive in its place if needed.
  private static Node? FixPragmaOnce(Parser.File file,
                                     bool needs_export_module) {
    // Collect the pragmas to fix.
    var pragma_onces = new List<PreprocessorDirective>();
    var blank_lines_after_pragma_once = new List<Text>();
    bool right_after_pragma_once = false;
    foreach (Node n in file.children) {
      if (n is PreprocessorDirective{ is_pragma_once: true } po) {
        pragma_onces.Add(po);
        right_after_pragma_once = true;
      } else if (right_after_pragma_once && n.is_blank_line) {
        blank_lines_after_pragma_once.Add((n as Text)!);
      } else {
        right_after_pragma_once = false;
      }
    }

    // And fix them.
    Node? node = null;
    foreach (PreprocessorDirective po in pragma_onces) {
      if (needs_export_module) {
        int position_in_parent = po.position_in_parent;
        node = new Export(parent: null);
        var module =
            new Module(parent: node, name: file.file_module_full_name);
        file.InsertChild(position_in_parent, node);
      } else {
        // If the `#pragma once` is not replaced by an `export module`, remove
        // the blank lines after it.
        foreach (Text blank in blank_lines_after_pragma_once) {
          file.RemoveChild(blank);
        }
      }
      file.RemoveChild(po);
    }
    return EnsureBlankLineAfter(node);
  }

  private static void FixProjectFile(FileInfo input_file,
                                     HashSet<string> extra_headers,
                                     bool dry_run) {
    string input_filename = input_file.FullName;
    string output_filename =
        input_file.DirectoryName + "\\" + input_file.Name + ".new";
    using (StreamReader reader = input_file.OpenText()) {
      using (StreamWriter writer = System.IO.File.CreateText(output_filename)) {
        while (!reader.EndOfStream) {
          string line = reader.ReadLine()!;
          bool is_extra_header = false;
          foreach (string extra_header in extra_headers) {
            if (line.Contains('"' +   extra_header + '"')) {
              // An extra header, it is now named .ixx.
              line = Regex.Replace(line, @"\.hpp"">", @".ixx"">");
              line = Regex.Replace(line, @"\.hpp"" />", @".ixx"" />");
              is_extra_header = true;
              break;
            }
          }
          if (!is_extra_header) {
            if (line.Contains("_body.hpp\">") ||
                line.Contains("_body.hpp\" />")) {
              // A body, don't change it.
            } else if (line.Contains(".hpp\">") || line.Contains(".hpp\" />")) {
              line = Regex.Replace(line, @"\.hpp"">", @".ixx"">");
              line = Regex.Replace(line, @"\.hpp"" />", @".ixx"" />");
            }
          }
          writer.WriteLine(line);
        }
      }
    }
    if (!dry_run) {
      System.IO.File.Move(output_filename, input_filename, overwrite: true);
    }
  }

  private static void SplitIncludes(List<Include> includes,
                                    HashSet<string> excluded,
                                    out SortedSet<Include> principia_includes,
                                    out SortedSet<Include> system_includes,
                                    out SortedSet<Include> third_party_includes,
                                    out SortedSet<Include>
                                        conditional_includes) {
    // Excluded headers are considered third-party headers, not Principia
    // headers.
    principia_includes =
        new SortedSet<Include>(includes.FindAll(inc =>
                                                    inc.is_principia &&
                                                    !inc.is_conditional &&
                                                    !excluded.Contains(
                                                        inc.path[^1] +
                                                        "." +
                                                        inc.extension)),
                               new IncludeComparer());
    system_includes =
        new SortedSet<Include>(
            includes.FindAll(inc => inc.is_system && !inc.is_conditional),
            new IncludeComparer());
    third_party_includes =
        new SortedSet<Include>(includes.FindAll(inc =>
                                                    (!inc.is_principia &&
                                                        !inc.is_system &&
                                                        !inc.is_conditional) ||
                                                    excluded.Contains(
                                                        inc.path[^1] +
                                                        "." +
                                                        inc.extension)),
                               new IncludeComparer());
    conditional_includes =
        new SortedSet<Include>(includes.FindAll(inc => inc.is_conditional),
                               new IncludeComparer());
  }

  public static void Run(string[] args) {
    // Parse the arguments.
    var projects = new List<DirectoryInfo>();
    var excluded = new HashSet<string>();
    var extra_headers = new HashSet<string>();
    bool dry_run = true;
    foreach (string arg in args) {
      if (arg.StartsWith("--") && arg.Contains(":")) {
        string[] split =
            arg.Split([ "--", ":" ], StringSplitOptions.None);
        string option = split[1];
        string value = split[2];
        if (option == "project") {
          projects.Add(new DirectoryInfo(value));
        } else if (option == "dry_run") {
          dry_run = bool.Parse(value);
        } else if (option == "exclude") {
          excluded.Add(value);
        } else if (option == "extra_header") {
          extra_headers.Add(value);
        } else {
          throw new ArgumentException("Unknown option " + option);
        }
      }
    }
    foreach (DirectoryInfo project in projects) {
      ProcessProject(project, excluded, extra_headers, dry_run);
    }
  }
}

} // namespace sourcerer
} // namespace principia
