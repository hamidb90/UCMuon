"""
ui_state.py — keeps persisted GUI settings alive while their widgets are hidden.

Streamlit drops the session state of a keyed widget that is not rendered in a
run (after one more rerun), so a control hidden by Basic mode, or by any
condition (the PARMA site while another spectrum is chosen), would lose its
value and save_settings() would then drop it from the autosave. This module
keeps a shadow copy of every persisted setting:

  restore_missing()  before any widget renders: put back the persisted keys
                     Streamlit dropped since the last run, and re-assign the
                     rest so a widget shown again gets its value in the browser;
  remember(values)   in save_settings(): update the copy with the current
                     values and return the full set to write;
  forget(*keys)      instead of st.session_state.pop() for a persisted key,
                     so the copy does not bring the old value back;
  forget_all()       for "Reset autosave".

The copy lives in st.session_state under _STORE (an ordinary key, never a
widget, so Streamlit never drops it).
"""
import streamlit as st

_STORE = "_ucmuon_persist_store"


def _store():
    return st.session_state.setdefault(_STORE, {})


def restore_missing():
    """Restore persisted keys that are missing from the session, and re-assign
    the others. Call before any widget renders (writing a widget key after its
    widget exists raises).

    The re-assignment matters in the browser: a widget that was hidden in the
    last runs (the Basic / Advanced switch) and is shown again would display
    its default, and send that back on the next interaction, overwriting the
    kept value (seen with Basic's detector checkbox, 2026-10-06). A value
    written through st.session_state is pushed to the browser."""
    ss = st.session_state
    for k, v in _store().items():
        ss[k] = ss[k] if k in ss else v


def remember(values):
    """Merge the current values into the copy and return all persisted values
    (current ones plus those of hidden widgets)."""
    store = _store()
    store.update(values)
    return dict(store)


def forget(*keys):
    """Drop persisted keys from the session and from the copy."""
    store = _store()
    for k in keys:
        st.session_state.pop(k, None)
        store.pop(k, None)


def forget_all():
    """Empty the copy (the session keys are the caller's business)."""
    st.session_state[_STORE] = {}
