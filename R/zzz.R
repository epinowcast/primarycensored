# Name attributes for base-R distribution wrappers used as primary event
# CDFs. The attr("name") convention lets .format_class() and
# .lookup_pprimary() identify functions that are not package-defined and
# therefore cannot be introspected via .extract_function_name().
#
# Note: stats::dunif and stats::punif already expose a ".Call(C_dunif, ...)"
# body, so .extract_function_name() returns the correct names for them
# without needing explicit attributes. No top-level assignments are needed
# for those functions.
