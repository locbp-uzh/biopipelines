# Colab usage

The LLM workflow talks to Google Colab through the official **Colab MCP server** (`github.com/googlecolab/colab-mcp`): the assistant creates notebooks, adds cells, and runs them on a Colab runtime via MCP tool calls. Do not use `colab-ssh`, `cloudflared`, `ngrok`, or any reverse-tunnel trick: Colab's terms prohibit remote shells on managed runtimes, and accounts have been suspended over it.

## One-time setup (Claude Code)

1. **Install `uv`** so `uv tool install` is available:
   ```
   python -m pip install uv
   ```
   Verify `uv --version` resolves.

2. **Install the server as a `uv` tool** — a self-contained venv with a pre-built `colab-mcp` executable:
   ```
   uv tool install git+https://github.com/googlecolab/colab-mcp
   ```
   On Windows this produces `%APPDATA%\uv\tools\colab-mcp\Scripts\colab-mcp.exe`. Use that path in step 4 — don't register the server as `uvx git+…`.

3. **Create the log directory the server writes to.** With `-l <dir>` the server crashes at startup (`FileNotFoundError` in `init_logger`) if `<dir>` doesn't exist — it does not create it. On Windows:
   ```
   mkdir %USERPROFILE%\colab-mcp-logs
   ```

4. **Register the MCP server at user scope** so it's available in every project and not committed to this repo. Point at the installed exe, not `uvx`:
   ```
   claude mcp add --scope user colab-mcp -- "<path-to>\colab-mcp.exe" -l "<path-to>\colab-mcp-logs"
   ```

   Verify with `claude mcp list` — it should show `colab-mcp … ✓ Connected`.

5. **Start a fresh Claude Code session.** MCP tools are injected into a session's tool list only at startup, so `colab-mcp` tools appear only in a conversation started *after* registration — relaunch-and-resume does not pick them up.

6. **Authorize with Google on first connect.** The server pulls in OAuth machinery (`fastmcp` auth + authlib JWT); expect a browser handshake the first time the agent calls `open_colab_browser_connection`. Sign in with the Google account whose Colab you want to use.

## Using Colab

The assistant calls `mcp__colab-mcp__open_colab_browser_connection` once per session to attach to a runtime. That call **unlocks the rest of the notebook-editing tools** — they're advertised post-connect rather than at session start, so the full tool surface only becomes visible after a successful connect. Confirm tool names against the live list as you use them.

The user picks the *runtime type* in the Colab UI; the assistant cannot request a GPU class.

What the MCP path does **not** give you: `userdata.get()` / Colab Secrets and `drive.mount()` are reportedly unavailable through it (true of the VS Code extension; assume the same here until shown otherwise). If a pipeline or `.install()` cell expects Drive, that's a gap.

## Workflow — running a test pipeline on Colab

1. Push the branch under test.
2. Via `colab-mcp`: create a notebook and add a first cell that clones the branch. The exact line depends on where the repo lives:

   **Public GitHub:**
   ```python
   !git clone -b <branch> https://github.com/<org>/biopipelines.git
   %cd biopipelines
   ```

   **Private repo over HTTPS** — use a token. In GitLab, the user has to create a *project access token* with role *`Reporter`* and scope *`read_repository`*. Both token *name* and *value* go into the clone URL:
   ```python
   from getpass import getpass
   tok_name = input("Token name: ")
   tok = getpass("Token value: ")
   !git clone -b <branch> https://{tok_name}:{tok}@gitlab.example.org/<group>/<project>.git
   %cd <project>
   ```
   `getpass` keeps the value out of the notebook output. Revoke the token when you're done and don't paste it into chat with the assistant either, since transcripts and MCP tool args get logged.

   Follow with the canonical Colab setup from `docs/user_manual.md` → "Google Colab".
3. Add a cell with the affected `.install()` call(s).
4. Add a minimal pipeline cell.
5. Execute cells in order; read cell outputs (success markers, traceback, `_log` files) directly via the MCP tools — no copy-paste round-trip with the user.
6. Iterate.

## Notes

- The settings in config.colab.yaml are automatically picked up.
- A Colab VM recycles at ~12 h (less on free, with idle disconnects). Keep test pipelines minimal.
- Guard Colab-specific fixes/patches that might break on cluster/local setups behind the scheduler (`if scheduler == "colab":`).
