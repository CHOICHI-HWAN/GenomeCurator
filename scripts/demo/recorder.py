"""Tiny screen recorder for the README GIFs: drives a headless Chrome over the DevTools protocol, shows a cursor and a
caption inside the page, grabs frames and assembles an optimised GIF. Needs `pip install -e ".[demo]"` and Google Chrome."""
import io
import json
import os
import shutil
import subprocess
import tempfile
import time
import urllib.request

import websocket
from PIL import Image, ImageChops

OVERLAY_JS = r"""
window.confirm = () => true; window.alert = () => {};
document.addEventListener('DOMContentLoaded', () => {
  const st = document.createElement('style');
  st.textContent = `#__cap{position:fixed;left:0;right:0;bottom:0;z-index:2147483000;background:rgba(20,30,45,.93);color:#fff;
   font:600 18px/1.4 'Noto Sans CJK KR','Noto Sans CJK HK',sans-serif;padding:12px 22px;display:none}
   #__cap b{color:#8fd0ff}
   #__cur{position:fixed;width:24px;height:24px;margin:-12px 0 0 -12px;border-radius:50%;background:rgba(255,80,60,.5);border:2px solid #fff;
   box-shadow:0 0 0 2px rgba(255,80,60,.95);z-index:2147483001;pointer-events:none;left:-60px;top:-60px;transition:left .3s ease,top .3s ease,transform .15s}
   #__cur.click{transform:scale(1.9);background:rgba(255,80,60,.15)}`;
  document.head.appendChild(st);
  const cap = document.createElement('div'); cap.id = '__cap'; document.body.appendChild(cap);
  const cur = document.createElement('div'); cur.id = '__cur'; document.body.appendChild(cur);
  window.__cap = (t) => { cap.innerHTML = t || ''; cap.style.display = t ? 'block' : 'none'; };
  window.__cur = (x, y, instant) => { cur.style.transition = instant ? 'none' : ''; cur.style.left = x + 'px'; cur.style.top = y + 'px'; };
  window.__click = (on) => cur.classList.toggle('click', on);
});
"""

FIND_JS = r"""
((css, text, nth) => {
  let els = [...document.querySelectorAll(css || 'button, a, summary, label, input[type=submit], [onclick]')];
  if (text) { const t = text.toLowerCase(); els = els.filter(e => (e.innerText || e.value || '').replace(/\s+/g, ' ').toLowerCase().includes(t)); }
  els = els.filter(e => { const r = e.getBoundingClientRect(); return r.width > 0 && r.height > 0; });
  const e = els[nth || 0];
  if (!e) return null;
  e.scrollIntoView({block: 'center', inline: 'center'});
  const r = e.getBoundingClientRect();
  window.__el = e;
  return {x: r.left + r.width / 2, y: r.top + r.height / 2, w: r.width, h: r.height};
})
"""


class Chrome:
    def __init__(self, width=1120, height=700, port=None):
        if port is None:   # a free port, so a Chrome left over from an earlier run can never be attached by mistake
            import socket
            with socket.socket() as sk:
                sk.bind(("127.0.0.1", 0))
                port = sk.getsockname()[1]
        exe = shutil.which("google-chrome") or shutil.which("chromium") or shutil.which("chrome")
        if not exe:
            raise SystemExit("Google Chrome / Chromium is required to record the demos")
        self.w, self.h = width, height
        self.profile = tempfile.mkdtemp(prefix="gc-demo-")
        self.proc = subprocess.Popen([exe, "--headless=new", "--no-sandbox", "--disable-gpu", "--hide-scrollbars", f"--remote-debugging-port={port}", "--remote-allow-origins=*",
                                      f"--user-data-dir={self.profile}", f"--window-size={width},{height}", "--force-device-scale-factor=1",
                                      "--lang=ko-KR", "about:blank"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.ws = None
        try:
            self._attach(port, width, height)
        except BaseException:
            self.close()   # never leave a headless Chrome behind when the attach fails
            raise

    def _attach(self, port, width, height):
        url = None
        for _ in range(100):
            try:
                tabs = json.load(urllib.request.urlopen(f"http://127.0.0.1:{port}/json/list", timeout=1))
                url = next(t["webSocketDebuggerUrl"] for t in tabs if t["type"] == "page")
                break
            except Exception:
                time.sleep(0.2)
        if not url:
            raise SystemExit("could not attach to Chrome")
        self.ws = websocket.create_connection(url, max_size=None, timeout=60)
        self._id = 0
        self.frames: list[tuple[Image.Image, int]] = []
        self.cx = self.cy = -60
        self.call("Page.enable")
        self.call("Emulation.setDeviceMetricsOverride", width=width, height=height, deviceScaleFactor=1, mobile=False)
        self.call("Page.addScriptToEvaluateOnNewDocument", source=OVERLAY_JS)

    # --- protocol
    def call(self, method, **params):
        self._id += 1
        self.ws.send(json.dumps({"id": self._id, "method": method, "params": params}))
        while True:
            msg = json.loads(self.ws.recv())
            if msg.get("id") == self._id:
                if "error" in msg:
                    raise RuntimeError(f"{method}: {msg['error']}")
                return msg.get("result", {})

    def js(self, expr):
        r = self.call("Runtime.evaluate", expression=expr, returnByValue=True, awaitPromise=True)
        if "exceptionDetails" in r:
            raise RuntimeError(r["exceptionDetails"].get("exception", {}).get("description", "js error"))
        return r.get("result", {}).get("value")

    def goto(self, url, wait=0.6):
        self.call("Page.navigate", url=url)
        self.wait("document.readyState==='complete'", 30)
        time.sleep(wait)
        self.js(f"window.__cur&&window.__cur({self.cx},{self.cy},true)")

    def wait(self, cond, timeout=30, every=0.25):
        end = time.time() + timeout
        while time.time() < end:
            try:
                if self.js(f"!!({cond})"):
                    return True
            except Exception:
                pass
            time.sleep(every)
        raise TimeoutError(cond)

    def settle(self, s=0.5):
        time.sleep(s)

    # --- showing things
    def caption(self, html):
        self.js(f"window.__cap&&window.__cap({json.dumps(html, ensure_ascii=False)})")

    def snap(self, ms=700):
        png = self.call("Page.captureScreenshot", format="png")["data"]
        import base64
        self.frames.append((Image.open(io.BytesIO(base64.b64decode(png))).convert("RGB"), ms))

    def hold(self, ms=1400):
        self.snap(ms)

    # --- acting
    def locate(self, css=None, text=None, nth=0):
        pos = self.js(f"{FIND_JS}({json.dumps(css)},{json.dumps(text)},{nth})")
        if not pos:
            raise LookupError(f"element not found: css={css!r} text={text!r}")
        return pos

    def move(self, x, y, frames=1):
        self.js(f"window.__cur({x},{y})")
        self.cx, self.cy = x, y
        time.sleep(0.4)
        for _ in range(frames):
            self.snap(260)

    def click(self, css=None, text=None, nth=0, after=0.5, show=True):
        pos = self.locate(css, text, nth)
        time.sleep(0.15)
        pos = self.locate(css, text, nth)
        if show:
            self.move(pos["x"], pos["y"], frames=1)
            self.js("window.__click(true)")
            self.snap(220)
            self.js("window.__click(false)")
        self.js("window.__el.click()")
        time.sleep(after)

    def type(self, css, text, step=3, ms=170):
        self.js(f"(()=>{{const e=document.querySelector({json.dumps(css)});e.focus();e.value='';}})()")
        for i in range(step, len(text) + step, step):
            part = text[:i]
            self.js(f"(()=>{{const e=document.querySelector({json.dumps(css)});e.value={json.dumps(part, ensure_ascii=False)};e.dispatchEvent(new Event('input',{{bubbles:true}}));}})()")
            self.snap(ms)

    def drag(self, x0, y0, x1, y1, steps=5):
        self.move(x0, y0)
        self.call("Input.dispatchMouseEvent", type="mouseMoved", x=x0, y=y0)
        self.call("Input.dispatchMouseEvent", type="mousePressed", x=x0, y=y0, button="left", buttons=1, clickCount=1)
        for i in range(1, steps + 1):
            x, y = x0 + (x1 - x0) * i / steps, y0 + (y1 - y0) * i / steps
            self.js(f"window.__cur({x},{y},true)")
            self.call("Input.dispatchMouseEvent", type="mouseMoved", x=x, y=y, button="left", buttons=1)
            time.sleep(0.12)
            self.snap(160)
        self.call("Input.dispatchMouseEvent", type="mouseReleased", x=x1, y=y1, button="left", buttons=0, clickCount=1)
        self.cx, self.cy = x1, y1

    def record_until(self, cond, every=1.0, timeout=180, ms=500, max_frames=24):
        """Grab frames while something runs (e.g. a progress modal); thin them out so the GIF stays small."""
        end, got = time.time() + timeout, []
        while time.time() < end:
            try:
                if self.js(f"!!({cond})"):
                    break
            except Exception:
                break          # page navigated away (the job redirected): that is the end
            self.snap(ms)
            got.append(len(self.frames) - 1)
            time.sleep(every)
        else:
            raise TimeoutError(cond)
        if len(got) > max_frames:   # keep an even subset
            keep = {got[int(i * (len(got) - 1) / (max_frames - 1))] for i in range(max_frames)}
            self.frames = [f for i, f in enumerate(self.frames) if i not in set(got) or i in keep]

    # --- output
    def save_gif(self, path, width=880, colors=112):
        out, durs, prev = [], [], None
        for im, ms in self.frames:
            h = round(im.height * width / im.width)
            sm = im.resize((width, h), Image.LANCZOS)
            if prev is not None and ImageChops.difference(prev, sm).getbbox() is None:
                durs[-1] += ms
                continue
            out.append(sm.quantize(colors=colors, method=Image.MEDIANCUT, dither=Image.NONE))
            durs.append(ms)
            prev = sm
        out[0].save(path, save_all=True, append_images=out[1:], duration=durs, loop=0, optimize=True, disposal=1)
        self.frames = []
        return len(out), os.path.getsize(path)

    def close(self):
        try:
            if self.ws:
                self.ws.close()
        finally:
            self.proc.terminate()
            shutil.rmtree(self.profile, ignore_errors=True)
