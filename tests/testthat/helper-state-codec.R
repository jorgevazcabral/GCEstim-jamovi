# Match jmvcore Element's state codec without the GUI-only protobuf bindings.
state_message <- function(state) {
  conn <- rawConnection(raw(), "r+")
  on.exit(close(conn))
  saveRDS(state, conn)
  list(state = memCompress(rawConnectionValue(conn)))
}
